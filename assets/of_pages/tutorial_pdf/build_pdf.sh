#!/usr/bin/env bash
# Build the tutorial PDF from OrthoFinder3.2_tutorial.html (and images/) with headless Chrome.
#   ./build_pdf.sh [-v VERSION] [output.pdf]
# -v: the release the PDF is for (e.g. 3.3): the version the text names (its
# <meta name="orthofinder-version">) is replaced by it in the install commands,
# download links and cover. Without -v, the text as it is.
# The tree figures in images/ are drawn from a run on ExampleData by draw_trees.py.
# Run by .github/workflows/update-release-pages.yml for each release.
set -euo pipefail
cd "$(dirname "$0")"
SRC=OrthoFinder3.2_tutorial.html
VERSION=""
if [[ "${1:-}" == "-v" ]]; then VERSION="${2:?-v needs a version}"; shift 2; fi
VERSION="${VERSION#v}"
OUT="${1:-OrthoFinder3.2_tutorial.pdf}"
[[ "$OUT" == /* ]] || OUT="$PWD/$OUT"
CHROME="$(command -v google-chrome || command -v google-chrome-stable || command -v chromium || command -v chromium-browser)" \
    || { echo "ERROR: Chrome or Chromium not found" >&2; exit 1; }

PAGE="$SRC"
if [[ -n "$VERSION" ]]; then
    # A copy next to the original (for the relative paths of images/), removed afterwards
    PAGE=".build_tutorial.html"
    trap 'rm -f "$PAGE"' EXIT
    python3 - "$SRC" "$PAGE" "$VERSION" <<'EOF'
import re, sys
src, out, new = sys.argv[1:]
text = open(src, encoding="utf-8").read()
old = re.search(r'<meta name="orthofinder-version" content="([^"]+)">', text).group(1)
if not re.fullmatch(r"\d+\.\d+(\.\d+)?", new):
    sys.exit("ERROR: not a release version: %s" % new)
o = re.escape(old)
# Only where the text names the release: v3.2, version 3.2, orthofinder=3.2,
# orthofinder-3.2.tar.gz, orthofinder-linux-intel-3.2.tar.gz, 3.2.post1.dev4
patterns = [r"(?<=v)%s(?![\d.]*\d)" % o, r"(?<=version )%s(?![\d.]*\d)" % o,
            r"(?<=orthofinder=)%s(?![\d.]*\d)" % o, r"(?<=-)%s(?=\.tar\.gz)" % o, r"%s(?=\.post1)" % o]
n = 0
for p in patterns:
    text, k = re.subn(p, new, text, flags=re.IGNORECASE)
    n += k
text = text.replace('content="%s"' % old, 'content="%s"' % new)
open(out, "w", encoding="utf-8").write(text)
print("Version %s -> %s in %d places" % (old, new, n))
EOF
fi

"$CHROME" --headless=new --disable-gpu --no-sandbox --no-pdf-header-footer \
    --run-all-compositor-stages-before-draw --virtual-time-budget=10000 \
    --generate-pdf-document-outline \
    --print-to-pdf="$OUT" "file://$PWD/$PAGE" 2>/dev/null
[[ -s "$OUT" ]] || { echo "ERROR: no PDF written" >&2; exit 1; }
echo "Wrote $OUT"
