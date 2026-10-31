#!/usr/bin/env python3
"""下载 Borkenhagen 2024 的 Data S1 序列和 Data S2 标签。"""

import argparse
from pathlib import Path
from urllib.request import Request, urlopen

DOI = "10.1111/irv.70044"
PMC_BASE = "https://pmc.ncbi.nlm.nih.gov/articles/instance/11634464/bin"
FILES = {
    "IRV-18-e70044-s003.txt": "DataS1_sequences.fasta",
    "IRV-18-e70044-s002.csv": "irv70044-sup-0003-supplementary_file2.csv",
}


def download(url: str, output: Path, timeout: int, overwrite: bool) -> None:
    if output.exists() and not overwrite:
        print(f"跳过已存在文件: {output}")
        return
    request = Request(url, headers={"User-Agent": "FluJump-release/0.1"})
    temporary = output.with_suffix(output.suffix + ".part")
    try:
        with urlopen(request, timeout=timeout) as response:
            data = response.read()
        if not data:
            raise RuntimeError(f"下载内容为空: {url}")
        temporary.write_bytes(data)
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"已下载: {output} ({output.stat().st_size:,} bytes)")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=f"获取 Borkenhagen 2024 补充数据（DOI: {DOI}）"
    )
    parser.add_argument("--output-dir", default="data/dataset_borkenhagen_raw")
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for remote_name, local_name in FILES.items():
        download(
            f"{PMC_BASE}/{remote_name}",
            output_dir / local_name,
            args.timeout,
            args.overwrite,
        )


if __name__ == "__main__":
    main()
