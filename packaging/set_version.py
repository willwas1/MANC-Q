"""Set MANC-Q's version number: python packaging/set_version.py 2.0.0

The number is kept in pyproject.toml only (MANC-Q reads it from there). This script changes it there and copies
it into the files that cannot read it themselves: CITATION.cff (with today's date), packaging/installer.iss and
the first line of constraints.txt. Without a number it only copies the current one."""
import datetime
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sub(rel, pattern, repl):
    p = os.path.join(ROOT, rel)
    with open(p, "rb") as fh:
        t = fh.read().decode("utf-8")
    t2, n = re.subn(pattern, repl, t, count=1, flags=re.M)
    if n != 1:
        sys.exit(f"{rel}: could not find the version to change ({pattern})")
    with open(p, "wb") as fh:
        fh.write(t2.encode("utf-8"))
    print(f"  {rel}")


def main(argv):
    with open(os.path.join(ROOT, "pyproject.toml"), encoding="utf-8") as fh:
        v = re.search(r'^version\s*=\s*"([^"]+)"', fh.read(), re.M).group(1)
    if argv:
        if not re.fullmatch(r"\d+\.\d+\.\d+", argv[0]):
            sys.exit(f"not a version number like 2.0.0: {argv[0]}")
        v = argv[0]
        sub("pyproject.toml", r'^version\s*=\s*"[^"]+"', f'version = "{v}"')
    sub("CITATION.cff", r"^version: .*$", f"version: {v}")
    sub("CITATION.cff", r"^date-released: .*$", f"date-released: {datetime.date.today():%Y-%m-%d}")
    sub("packaging/installer.iss", r'#define MyAppVersion "[^"]+"', f'#define MyAppVersion "{v}"')
    sub("constraints.txt", r"^# Exact package versions MANC-Q \S+ was tested",
        f"# Exact package versions MANC-Q {v} was tested")
    print(f"MANC-Q version {v}")


if __name__ == "__main__":
    main(sys.argv[1:])
