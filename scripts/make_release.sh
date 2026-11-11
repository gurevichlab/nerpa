#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

release_name="${1:-}"
git_ref="${2:-HEAD}"

if [[ ! "$release_name" =~ ^nerpa-ms-[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "Usage: $0 nerpa-ms-X.Y.Z [GIT_REF]" >&2
    exit 1
fi

if ! git rev-parse --verify --quiet "${git_ref}^{commit}" >/dev/null; then
    echo "Error: Git ref '$git_ref' does not identify a commit." >&2
    exit 1
fi

output_dir="$repo_root/dist"
archive_path="$output_dir/$release_name.tar.gz"
temporary_tar="$(mktemp)"

trap 'rm -f "$temporary_tar"' EXIT

mkdir -p "$output_dir"

# git archive includes only committed files from the selected revision.
release_commit="$(git rev-parse "${git_ref}^{commit}")"
git archive \
    --format=tar \
    --prefix="$release_name/" \
    "$release_commit" > "$temporary_tar"

# Suppress the gzip timestamp so rebuilding the same tag is reproducible.
gzip -n -c "$temporary_tar" > "$archive_path"
(
    cd "$output_dir"
    sha256sum "$release_name.tar.gz" > "$release_name.tar.gz.sha256"
)

echo "Archived commit: $release_commit"
echo "Created:"
echo "  $archive_path"
echo "  $archive_path.sha256"