#!/usr/bin/env python3
"""
Mirror a GitHub issue from a core repo into a consumer repo.

 - copies title, body, labels, assignees, milestone, issue type
 - marks the consumer issue as "blocked by" the core issue (native GitHub dependency)
 - adds it to the same Projects v2 board(s) and copies the project field values
 - idempotent: re-running updates the existing mirror instead of creating a new one

Requirements: gh CLI (recent), authenticated with the `project` scope:
    gh auth refresh -s project

Usage:
    ./mirror_issue.py 556                       # mirrors ottolink-dev/HighMap#556 -> ottolink-dev/Hesiod
    ./mirror_issue.py ORG/core#42 --to ORG/app
    ./mirror_issue.py https://github.com/ORG/core/issues/42 --to ORG/app --title-prefix "[core] "
    ./mirror_issue.py ORG/core#42 --to ORG/app --sync-state --skip-field Status
"""
import argparse
import json
import re
import subprocess
import sys

DEFAULT_CORE_REPO = "ottolink-dev/HighMap"
DEFAULT_CONSUMER_REPO = "ottolink-dev/Hesiod"
MARK = "core-mirror"    # hidden marker in the body


class GhError(Exception):
    pass


def run(args, payload=None):
    r = subprocess.run(
        ["gh", *args],
        input=None if payload is None else json.dumps(payload),
        capture_output=True, text=True,
    )
    if r.returncode:
        raise GhError((r.stderr or r.stdout).strip())
    return r.stdout


def api(method, path, payload=None):
    args = ["api", "-X", method, path]
    if payload is not None:
        args += ["--input", "-"]
    out = run(args, payload)
    return json.loads(out) if out.strip() else None


def gql(query, **variables):
    out = json.loads(run(["api", "graphql", "--input", "-"],
                         {"query": query, "variables": variables}))
    if out.get("errors"):
        raise GhError(json.dumps(out["errors"]))
    return out["data"]


def parse_ref(s):
    s = s.strip()
    if s.isdigit() or (s.startswith("#") and s[1:].isdigit()):
        num = int(s.lstrip("#"))
        o, r = DEFAULT_CORE_REPO.split("/", 1)
        return o, r, num

    m = re.search(r"(?:github\.com/)?([\w.-]+)/([\w.-]+)(?:/issues/|#)(\d+)", s)
    if not m:
        sys.exit(f"cannot parse issue reference: {s}")
    return m.group(1), m.group(2), int(m.group(3))


FETCH = """
query($o:String!,$r:String!,$n:Int!){
  repository(owner:$o,name:$r){
    issue(number:$n){
      id databaseId number title body state stateReason url
      author{login}
      labels(first:100){nodes{name color description}}
      assignees(first:20){nodes{login}}
      milestone{title description dueOn}
      issueType{name}
      projectItems(first:20){nodes{
        project{id number title}
        fieldValues(first:50){nodes{
          __typename
          ... on ProjectV2ItemFieldTextValue{text field{... on ProjectV2FieldCommon{id name}}}
          ... on ProjectV2ItemFieldNumberValue{number field{... on ProjectV2FieldCommon{id name}}}
          ... on ProjectV2ItemFieldDateValue{date field{... on ProjectV2FieldCommon{id name}}}
          ... on ProjectV2ItemFieldSingleSelectValue{optionId field{... on ProjectV2FieldCommon{id name}}}
          ... on ProjectV2ItemFieldIterationValue{iterationId field{... on ProjectV2FieldCommon{id name}}}
        }}
      }}
    }
  }
}
"""

ADD_ITEM = """
mutation($p:ID!,$c:ID!){
  addProjectV2ItemById(input:{projectId:$p,contentId:$c}){item{id}}
}
"""

SET_FIELD = """
mutation($p:ID!,$i:ID!,$f:ID!,$v:ProjectV2FieldValue!){
  updateProjectV2ItemFieldValue(input:{projectId:$p,itemId:$i,fieldId:$f,value:$v}){
    projectV2Item{id}
  }
}
"""


def to_value(fv):
    t = fv.get("__typename")
    if t == "ProjectV2ItemFieldTextValue":
        return {"text": fv["text"]}
    if t == "ProjectV2ItemFieldNumberValue":
        return {"number": fv["number"]}
    if t == "ProjectV2ItemFieldDateValue":
        return {"date": fv["date"]}
    if t == "ProjectV2ItemFieldSingleSelectValue":
        return {"singleSelectOptionId": fv["optionId"]}
    if t == "ProjectV2ItemFieldIterationValue":
        return {"iterationId": fv["iterationId"]}
    return None  # built-in fields (title, labels, assignees, ...) are not settable here


def ensure_labels(repo, labels):
    for l in labels:
        name = l.get("name") if isinstance(l, dict) else str(l)
        if not name:
            continue
        try:
            api("POST", f"repos/{repo}/labels", {
                "name": name, "color": l.get("color") if isinstance(l, dict) and l.get("color") else "ededed",
                "description": l.get("description") if isinstance(l, dict) and l.get("description") else "",
            })
        except GhError:
            pass  # already exists


def ensure_milestone(repo, ms):
    if not ms:
        return None
    pages = json.loads(run(["api", "--paginate", "--slurp",
                            f"repos/{repo}/milestones?state=all&per_page=100"]))
    for page in pages:
        for m in page:
            if m["title"] == ms["title"]:
                return m["number"]
    payload = {"title": ms["title"], "description": ms.get("description") or ""}
    if ms.get("dueOn"):
        payload["due_on"] = ms["dueOn"]
    return api("POST", f"repos/{repo}/milestones", payload)["number"]


def find_mirror(consumer, tag):
    out = run(["issue", "list", "--repo", consumer, "--state", "all",
               "--limit", "1000", "--json", "number,body,id"])
    for i in json.loads(out):
        if tag in (i.get("body") or ""):
            return i
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("issue", help="core issue: issue number (assumes HighMap), ORG/repo#N, or full URL")
    ap.add_argument("--to", default=DEFAULT_CONSUMER_REPO,
                    help=f"consumer repo: ORG/repo (default: {DEFAULT_CONSUMER_REPO})")
    ap.add_argument("--title-prefix", default="")
    ap.add_argument("--project", type=int, help="only sync this project number (default: all shared)")
    ap.add_argument("--skip-field", action="append", default=[], help="project field name to not copy (repeatable)")
    ap.add_argument("--sync-state", action="store_true", help="also mirror open/closed state")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    o, r, n = parse_ref(args.issue)
    core_ref = f"{o}/{r}#{n}"
    consumer = args.to

    core = gql(FETCH, o=o, r=r, n=n)["repository"]["issue"]
    if not core:
        sys.exit(f"issue not found: {core_ref}")

    labels = core["labels"]["nodes"]
    label_names = [l["name"] for l in labels if l.get("name")]
    assignees = [a["login"] for a in core["assignees"]["nodes"]]
    tag = f"<!-- {MARK}: {core_ref} -->"
    author = (core.get("author") or {}).get("login", "unknown")
    body = (
        f"{tag}\n"
        f"> Mirror of {core['url']} (opened by {author}).\n"
        f"> Blocked by {core_ref}. Edit the core issue, not this one.\n\n"
        f"{core.get('body') or ''}"
    )
    title = args.title_prefix + core["title"]

    if args.dry_run:
        print(json.dumps({"title": title, "labels": label_names, "assignees": assignees,
                          "milestone": core["milestone"], "type": core["issueType"],
                          "projects": [i["project"]["title"] for i in core["projectItems"]["nodes"]],
                          "body": body}, indent=2))
        return

    ensure_labels(consumer, labels)
    ms_number = ensure_milestone(consumer, core["milestone"])

    existing = find_mirror(consumer, tag)
    if existing:
        num = existing["number"]
        patch = {"title": title, "body": body, "assignees": assignees, "milestone": ms_number}
        if args.sync_state:
            if core["state"] == "CLOSED":
                patch["state"] = "closed"
                patch["state_reason"] = (core.get("stateReason") or "completed").lower()
            else:
                patch["state"] = "open"
        issue = api("PATCH", f"repos/{consumer}/issues/{num}", patch)
        api("POST", f"repos/{consumer}/issues/{num}/labels", {"labels": label_names})
        print(f"updated {consumer}#{num}")
    else:
        payload = {"title": title, "body": body, "labels": label_names, "assignees": assignees}
        if ms_number:
            payload["milestone"] = ms_number
        if core["issueType"]:
            payload["type"] = core["issueType"]["name"]
        try:
            issue = api("POST", f"repos/{consumer}/issues", payload)
        except GhError as e:
            if "type" not in payload:
                raise
            print(f"warning: issue type not applied ({e}); retrying without it", file=sys.stderr)
            payload.pop("type")
            issue = api("POST", f"repos/{consumer}/issues", payload)
        num = issue["number"]
        print(f"created {issue['html_url']}")

    # consumer issue is "blocked by" the core issue (native dependency)
    try:
        api("POST", f"repos/{consumer}/issues/{num}/dependencies/blocked_by",
            {"issue_id": core["databaseId"]})
        print(f"{consumer}#{num} is now blocked by {core_ref}")
    except GhError as e:
        print(f"note: blocked-by not added ({e}); the body still carries the reference", file=sys.stderr)

    # project membership + field values (same project => same field/option/iteration ids)
    node_id = issue["node_id"]
    for item in core["projectItems"]["nodes"]:
        proj = item["project"]
        if args.project and proj["number"] != args.project:
            continue
        new_item = gql(ADD_ITEM, p=proj["id"], c=node_id)["addProjectV2ItemById"]["item"]["id"]
        for fv in item["fieldValues"]["nodes"]:
            value = to_value(fv)
            field = fv.get("field") or {}
            if value is None or not field or field["name"] in args.skip_field:
                continue
            try:
                gql(SET_FIELD, p=proj["id"], i=new_item, f=field["id"], v=value)
                print(f"  [{proj['title']}] {field['name']} set")
            except GhError as e:
                print(f"  [{proj['title']}] {field['name']} skipped: {e}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except GhError as e:
        sys.exit(f"gh error: {e}")
