#!/bin/bash
# build the independent proof checker drat-trim (github.com/marijnheule/drat-trim) into tools/
set -euo pipefail
cd "$(dirname "$0")/.."; mkdir -p tools
[ -d tools/drat-trim ] || git clone -q https://github.com/marijnheule/drat-trim.git tools/drat-trim
make -s -C tools/drat-trim >/dev/null 2>&1 || make -C tools/drat-trim
ls -la tools/drat-trim/drat-trim && echo "drat-trim ready"
