#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

version_file="${1:-}"

if [[ -z "$version_file" ]]; then
    echo "Usage: $0 VERSION_FILE" >&2
    exit 1
fi

release_commit="$(git rev-parse HEAD)"
short_commit="$(git rev-parse --short=12 HEAD)"
version="$(tr -d '[:space:]' < "$version_file")"

if [[ ! "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "Error: '$version_file' does not contain an X.Y.Z version." >&2
    exit 1
fi

release_tag="nerpa-ms-$version"
tag_commit="$(git rev-parse --verify --quiet "refs/tags/$release_tag^{commit}" || true)"

if [[ "$tag_commit" == "$release_commit" ]]; then
    release_name="$release_tag"
else
    release_name="$release_tag-$short_commit"
fi

output_dir="$repo_root/dist"
archive_path="$output_dir/$release_name.tar.gz"
temporary_tar="$(mktemp)"

trap 'rm -f "$temporary_tar"' EXIT

mkdir -p "$output_dir"

# Archive committed files only; .gitattributes controls release exclusions.
git archive \
    --format=tar \
    --prefix="$release_name/" \
    "$release_commit" > "$temporary_tar"

# Suppress the gzip timestamp so rebuilding the same commit is reproducible.
gzip -n -c "$temporary_tar" > "$archive_path"
(
    cd "$output_dir"
    sha256sum "$release_name.tar.gz" > "$release_name.tar.gz.sha256"
)

echo "Version: $version"
echo "Archived commit: $release_commit"

if [[ "$release_name" == "$release_tag" ]]; then
    echo "Release tag: $release_tag"
else
    echo "Release tag '$release_tag' is not attached to HEAD; created a development archive."
fi

echo "Created:"
echo "  $archive_path"
echo "  $archive_path.sha256"