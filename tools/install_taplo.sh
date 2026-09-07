#!/usr/bin/env bash
# Pinned Linux x86_64 release used by local checks and CI.
set -euo pipefail
install_dir="${1:?Usage: install_taplo.sh INSTALL_DIRECTORY}"
mkdir -p "$install_dir"
archive="$install_dir/taplo-0.10.0.gz"
curl --fail --silent --show-error --location \
  https://github.com/tamasfe/taplo/releases/download/0.10.0/taplo-linux-x86_64.gz \
  --output "$archive"
expected=8fe196b894ccf9072f98d4e1013a180306e17d244830b03986ee5e8eabeb6156
actual=$(sha256sum "$archive")
[[ "${actual%% *}" == "$expected" ]] || { echo 'Taplo checksum mismatch' >&2; exit 1; }
gzip -dc "$archive" > "$install_dir/taplo"
chmod +x "$install_dir/taplo"
"$install_dir/taplo" --version
