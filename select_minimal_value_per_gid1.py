#!/usr/bin/env python
# usage: python select_gid_pair_with_min_per_gid.py input.tsv output.tsv

import sys

def is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False

def main(infile, outfile):
    header = None
    records = []  # (value, idx, key1, key2, line)
    seen_key1 = set()

    # 读取并准备数据
    with open(infile, "r", encoding="utf-8") as f:
        for idx, raw in enumerate(f):
            line = raw.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue

            # 自动识别表头（第三列不是数值则认为是表头）
            if idx == 0 and not is_float(cols[2]):
                header = line
                continue

            k1, k2 = cols[0], cols[1]
            vstr = cols[2]
            if not is_float(vstr):
                continue
            v = float(vstr)

            records.append((v, idx, k1, k2, line))
            seen_key1.add(k1)

    # 按第3列值升序（并以输入顺序 idx 保持稳定）
    records.sort(key=lambda t: (t[0], t[1]))

    used_key1 = set()
    used_key2 = set()
    kept = []

    # 贪心选择：从最小值开始，只要 key1 和 key2 都未使用就保留
    for v, idx, k1, k2, line in records:
        if k1 in used_key1 or k2 in used_key2:
            continue
        kept.append(line)
        used_key1.add(k1)
        used_key2.add(k2)

    # 写输出文件（保持已选行的升序）
    with open(outfile, "w", encoding="utf-8") as out:
        if header is not None:
            out.write(header + "\n")
        for line in kept:
            out.write(line + "\n")

    # 未被使用的 key1 打到标准错误
    unused_key1 = sorted(seen_key1 - used_key1)
    if unused_key1:
        sys.stderr.write(f"Unused key1 count: {len(unused_key1)}\n")
        for k in unused_key1:
            sys.stderr.write(str(k) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python greedy_unique_min.py input.tsv output.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
