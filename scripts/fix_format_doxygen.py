#!/usr/bin/env python3
"""Clean up Doxygen comments.

1. Align continuation lines of @param/@return/... with the description column.
2. Pull wrapped '// ...' continuations back onto their '///<' line.
3. Re-join statements that were split across lines by a trailing '///<'.
"""
import argparse
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------- block comments
# @param[in] name   description  /  @tparam T  desc  /  @retval 0  desc
NAMED = re.compile(
    r'^(\s*\*\s*)(@(?:param(?:\[[^\]]*\])?|tparam|retval|throws|exception)\s+\S+\s+)(\S.*)?$')
# @return   description
UNNAMED = re.compile(r'^(\s*\*\s*)(@(?:returns?|result)\s+)(\S.*)?$')
ANY_TAG = re.compile(r'^\s*\*\s*[@\\]\w+')
BLANK = re.compile(r'^\s*\*\s*$')
STAR = re.compile(r'^(\s*\*)\s*(.*)$')

# ---------------------------------------------------------------- trailing comments
# code ... ///< text   (also //!<)
TRAIL = re.compile(r'^(.*\S\s*(?:///|//!)<\s*)(\S.*?)\s*$')
# plain "// text" line (not /// or //!)
CONT = re.compile(r'^\s*//(?![/!])\s?(.*\S)\s*$')
SENTENCE_END = ('.', '!', '?')

# statement ending with a ///< comment
TRAILING = re.compile(r'^(\s*\S.*?[;,}])\s*///<\s*(\S.*?)\s*$')
OPENERS, CLOSERS = '({[', ')}]'
ALIGN = re.compile(r'^(\s*[^\s/].*?)\s*///<\s*(.*?)\s*$')

# ================================================================ block comments
def tag_match(line):
    return NAMED.match(line) or UNNAMED.match(line)


def fix_block(block):
    # Description column = widest "prefix + tag + name" among tag lines
    cols = [len(m.group(1)) + len(m.group(2))
            for l in block if (m := tag_match(l))]
    if not cols:
        return block
    col = max(cols)

    out, in_tag = [], False
    for line in block:
        if tag_match(line):
            in_tag = True
            out.append(line)
        elif BLANK.match(line) or ANY_TAG.match(line) or '*/' in line:
            in_tag = False
            out.append(line)
        elif in_tag and (m := STAR.match(line)):
            star, text = m.groups()
            out.append(star + ' ' * max(1, col - len(star)) + text)
        else:
            out.append(line)
    return out


# ================================================================ trailing comments
def join_trailing_comments(text):
    """Pull wrapped '// ...' continuations back onto their '///<' line."""
    lines = text.split('\n')
    out, i = [], 0
    while i < len(lines):
        m = TRAIL.match(lines[i])
        if not m:
            out.append(lines[i])
            i += 1
            continue
        head, body = m.groups()
        j = i + 1
        while (not body.endswith(SENTENCE_END) and j < len(lines)
               and (c := CONT.match(lines[j]))):
            body += ' ' + c.group(1)
            j += 1
        out.append(head + body)
        i = j
    return '\n'.join(out)


def align_trailing(text, gap=1, max_col=100):
    """Align the '///<' of consecutive lines to one column."""
    lines = text.split('\n')
    ms = [ALIGN.match(l) for l in lines]
    out, i = [], 0
    while i < len(lines):
        if not ms[i]:
            out.append(lines[i])
            i += 1
            continue
        j = i
        while j < len(lines) and ms[j]:
            j += 1
        ok = {k for k in range(i, j) if len(ms[k].group(1)) + gap <= max_col}
        col = max((len(ms[k].group(1)) for k in ok), default=0) + gap
        for k in range(i, j):
            if k in ok:
                code, comment = ms[k].groups()
                out.append(code.ljust(col) + '///< ' + comment)
            else:
                out.append(lines[k])
        i = j
    return '\n'.join(out)


def _code_only(line):
    return re.sub(r'\s*//.*$', '', line)


def _is_boundary(line):
    s = _code_only(line).strip()      # ignore trailing comments
    return (not s                      # blank or comment-only line
            or s.startswith(('#', '/*', '*'))
            or s.endswith((';', '{', '}', ':', '*/')))


def _stmt_start(lines, end, max_lines=20):
    """Walk back from `end` to the first line of the statement."""
    bal = 0
    for k in range(end, max(end - max_lines, -1), -1):
        code = _code_only(lines[k])
        bal += sum(map(code.count, OPENERS)) - sum(map(code.count, CLOSERS))
        if bal == 0 and (k == 0 or _is_boundary(lines[k - 1])):
            return k
    return None


def fix_trailing(text):
    """Re-join statements broken across lines by a trailing ///< comment."""
    lines = text.split('\n')
    out, i = [], 0
    while i < len(lines):
        m = TRAILING.match(lines[i])
        start = _stmt_start(lines, i) if m else None
        if start is None or start == i:              # no match, or already one line
            out.append(lines[i])
            i += 1
            continue
        if any('///<' in l for l in lines[start:i]):  # never merge other members
            out.append(lines[i])
            i += 1
            continue
        del out[len(out) - (i - start):]             # drop raw lines already emitted
        code, comment = m.groups()
        parts = [lines[start].rstrip()] + [l.strip() for l in lines[start + 1:i]]
        parts.append(code.strip())
        out.append(f'{" ".join(parts)} ///< {comment}')
        i += 1
    return '\n'.join(out)


# ================================================================ driver
def process(text):
    lines = text.split('\n')
    out, block, in_doc = [], [], False
    for line in lines:
        if not in_doc:
            out.append(line)
            if '/**' in line and '*/' not in line:
                in_doc, block = True, []
            continue
        block.append(line)
        if '*/' in line:
            out.extend(fix_block(block))
            in_doc = False
    if in_doc:  # unterminated comment: leave untouched
        out.extend(block)
    return align_trailing(fix_trailing(join_trailing_comments('\n'.join(out))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+', type=Path)
    ap.add_argument('--check', action='store_true',
                    help='report only, exit 1 if changes needed')
    args = ap.parse_args()

    dirty = False
    for f in args.files:
        src = f.read_text(encoding='utf-8')
        new = process(src)
        if new != src:
            dirty = True
            print(f'{"would fix" if args.check else "fixed"}: {f}')
            if not args.check:
                f.write_text(new, encoding='utf-8')
    sys.exit(1 if (dirty and args.check) else 0)


if __name__ == '__main__':
    main()
