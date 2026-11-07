#!/usr/bin/env python3
"""Format or check production C++ sources under include/ and src/.

Use --clang-format or CLANG_FORMAT to select a clang-format executable.
Use the same formatter version across contributors for identical wrapping;
the current style was checked with clang-format 23.1.0.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess


def source_files(root, scope):
    if scope == "solver":
        files = [root / "include/MonteCarloSolver.h", root / "src/MonteCarloSolver.cpp"]
        folders = [root / "include/solver", root / "src/solver"]
    else:
        files = []
        folders = [root / "include", root / "src"]
    for folder in folders:
        files.extend(path for path in folder.rglob("*") if path.suffix in {".h", ".cpp"})
    return sorted(set(files))


def main(default_scope="all"):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Report formatting differences only")
    parser.add_argument("--scope", choices=("all", "solver"), default=default_scope)
    parser.add_argument("--clang-format", default=os.environ.get("CLANG_FORMAT", "clang-format"))
    args = parser.parse_args()
    executable = shutil.which(args.clang_format)
    if executable is None:
        parser.error("clang-format not found; set CLANG_FORMAT or pass --clang-format")

    root = Path(__file__).resolve().parents[1]
    files = source_files(root, args.scope)
    flags = ["--dry-run", "--Werror"] if args.check else ["-i"]
    result = subprocess.run([executable, "--style=file", *flags, *map(str, files)], cwd=root)
    if result.returncode == 0:
        print(f"{'Checked' if args.check else 'Formatted'} {len(files)} C++ files ({args.scope})")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
