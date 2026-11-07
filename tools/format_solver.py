#!/usr/bin/env python3
"""Compatibility entry point for formatting only the solver layer."""

from format_cpp import main


if __name__ == "__main__":
    raise SystemExit(main(default_scope="solver"))
