"""
Update the website (the gh-pages branch) for a new OrthoFinder release:
  - docs/download.markdown: the download button, its download count and the
    install commands name the new version
  - _posts/<date>-OrthoFinder-v<version>.markdown: the release's post, with its
    download button and release notes (the Older Versions page lists the posts)

Run by .github/workflows/update-release-pages.yml when a release is published.
Running it again for the same release changes nothing (an existing post is kept).

    python update_release_pages.py SITE_DIR VERSION [--date YYYY-MM-DD] [--notes-file FILE]
"""

import argparse
import datetime
import os
import re
import sys

REPO_URL = "https://github.com/OrthoFinder/OrthoFinder"
# A version in a release asset name or tag: 3.1.5, 3.2, 3.2.0rc1 ...
VERSION_RE = r"\d+\.\d+(?:\.\d+)?(?:[a-z]+\d*)?"

DOWNLOAD_ITEM = """<div class="download-item">
  📦
    <a href="{url}" class="btn">
        {asset}
    </a>
    <span class="download-count" data-asset="{asset}">
    <i class="fa fa-download" aria-hidden="true"></i>
    </span>
</div>"""

POST = """---
layout: post
title:  "OrthoFinder-v{version}"
date:   {date} 18:00:00 +0000
categories: release
permalink: /:categories/:day/:month/:year/:title.html
---

{download_item}
<!--more-->

{notes}
"""


def standalone_asset(version):
    return "orthofinder-linux-intel-%s.tar.gz" % version


def asset_url(version):
    return "%s/releases/download/v%s/%s" % (REPO_URL, version, standalone_asset(version))


def update_download_page(text, version):
    """The download page's text with every release it names replaced by version."""
    text = re.sub(r"/releases/download/v%s/" % VERSION_RE, "/releases/download/v%s/" % version, text)
    text = re.sub(r"orthofinder-linux-intel-%s\.tar\.gz" % VERSION_RE, standalone_asset(version), text)
    # pip install from GitHub: the release, not the development code on main
    text = re.sub(r"(pip install git\+%s\.git)(@v%s)?" % (re.escape(REPO_URL), VERSION_RE),
                  r"\1@v%s" % version, text)
    return text


def write_post(site_dir, version, date, notes):
    """Write the release's post unless it already has one. Returns its path, or None if kept."""
    posts_dir = os.path.join(site_dir, "_posts")
    os.makedirs(posts_dir, exist_ok=True)
    suffix = "-OrthoFinder-v%s.markdown" % version
    existing = [fn for fn in os.listdir(posts_dir) if fn.endswith(suffix)]
    if existing:
        return None
    fn = os.path.join(posts_dir, date + suffix)
    download_item = DOWNLOAD_ITEM.format(url=asset_url(version), asset=standalone_asset(version))
    # The release notes are Markdown (with CRLF if written in GitHub's web editor)
    notes = notes.replace("\r\n", "\n").strip() or "OrthoFinder v%s." % version
    text = POST.format(version=version, date=date, download_item=download_item, notes=notes)
    with open(fn, "w", encoding="utf-8", newline="\n") as outfile:
        outfile.write(text)
    return fn


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("site_dir", help="a checkout of the gh-pages branch")
    parser.add_argument("version", help="the release, e.g. 3.2 (a leading v is ignored)")
    parser.add_argument("--date", default=datetime.date.today().isoformat(), help="release date, YYYY-MM-DD")
    parser.add_argument("--notes-file", help="the release notes (Markdown); default: a one-line note")
    args = parser.parse_args(argv)

    version = args.version[1:] if args.version.startswith("v") else args.version
    if not re.fullmatch(VERSION_RE, version):
        sys.exit("ERROR: not a release version: %s" % args.version)
    datetime.date.fromisoformat(args.date)       # (fails on a malformed date)
    notes = ""
    if args.notes_file:
        with open(args.notes_file, encoding="utf-8") as infile:
            notes = infile.read()

    download_fn = os.path.join(args.site_dir, "docs", "download.markdown")
    with open(download_fn, encoding="utf-8") as infile:
        text = infile.read()
    new_text = update_download_page(text, version)
    if new_text != text:
        with open(download_fn, "w", encoding="utf-8", newline="\n") as outfile:
            outfile.write(new_text)
        print("Updated %s" % download_fn)
    else:
        print("Unchanged %s" % download_fn)

    post_fn = write_post(args.site_dir, version, args.date, notes)
    print("Wrote %s" % post_fn if post_fn else "Kept the existing post for v%s" % version)


if __name__ == "__main__":
    main()
