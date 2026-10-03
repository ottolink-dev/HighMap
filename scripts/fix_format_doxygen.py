#!/usr/bin/env python3
"""Align Doxygen continuation lines with the description column of their tag."""
import argparse
import re
import sys
from pathlib import Path

# @param[in] name   description  /  @tparam T  desc  /  @retval 0  desc
NAMED = re.compile(
    r'^(\s*\*\s*)(@(?:param(?:\[[^\]]*\])?|tparam|retval|throws|exception)\s+\S+\s+)(\S.*)?$')
# @return   description
UNNAMED = re.compile(r'^(\s*\*\s*)(@(?:returns?|result)\s+)(\S.*)?$')
ANY_TAG = re.compile(r'^\s*\*\s*[@\\]\w+')
BLANK = re.compile(r'^\s*\*\s*$')
STAR = re.compile(r'^(\s*\*)\s*(.*)$')

# code ... ///< text   (also //!<)
TRAIL = re.compile(r'^(.*\S\s*(?:///|//!)<\s*)(\S.*?)\s*$')
# plain "// text" line (not /// or //!)
CONT = re.compile(r'^\s*//(?![/!])\s?(.*\S)\s*$')
SENTENCE_END = ('.', '!', '?')


def tag_match(line):
    return NAMED.match(line) or UNNAMED.match(line)


def fix_block(block):
    # Description column of the block = widest "prefix + tag + name" among tag lines
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
    return join_trailing_comments('\n'.join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+', type=Path)
    ap.add_argument('--check', action='store_true', help='report only, exit 1 if changes needed')
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
