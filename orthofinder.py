#!/usr/bin/env python3
import sys
import os

orthofinder_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "src"))
sys.path.insert(0, orthofinder_dir)

if __name__ == "__main__":
    from orthofinder.run.__main__ import main
    args = sys.argv[1:]
    main(args)