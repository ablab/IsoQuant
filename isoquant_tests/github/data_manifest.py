#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Track and verify the integrity of CI test input data.

CI test inputs (genomes, BAMs, annotations, ground-truth tables) live outside the
repository on shared storage, so git has no record of them. This tool records each
input's size and SHA-256 in a manifest that *is* committed, which gives two things
the configs alone cannot:

* provenance - a baseline can be tied to the exact input bytes that produced it;
* integrity  - a file replaced or truncated under CI is caught immediately, rather
               than showing up later as an unexplained baseline drift.

Usage:
    data_manifest.py generate [--configs DIR] [--manifest FILE]
    data_manifest.py verify   [--configs DIR] [--manifest FILE] [--deep]
    data_manifest.py verify   --config-file path/to/One.yaml [--deep]

`verify` checks existence and size by default, which costs a stat() per file.
`--deep` also re-hashes, which reads every byte and is meant for a scheduled job
rather than for every test run.
"""

import argparse
import glob
import hashlib
import logging
import os
import sys
from typing import Dict, Iterable, List, Optional, Set, Tuple

import yaml

# Run by absolute path from CI, so the repo root is not on sys.path by default;
# prepend it so isoquant_lib resolves to this checkout, not an installed copy.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..")))

from isoquant_lib.utils.error_codes import IsoQuantExitCode

log = logging.getLogger('DataManifest')

# Config keys holding a single path.
SINGLE_PATH_KEYS: Set[str] = {
    "genome", "genedb", "reduced_db", "reference_tpm", "reference_gene_tpm",
    "yaml", "molecule", "reference_polya_gtf", "truth_file",
}
# Config keys holding one or more space-separated paths.
MULTI_PATH_KEYS: Set[str] = {"bam", "ubam", "reads", "barcodes"}
# Keys pointing at directories or outputs, deliberately excluded: `resume`
# checkpoints are regenerated, and `output` is written, not read.
EXCLUDED_KEYS: Set[str] = {"resume", "output"}

# `reduced_db` and similar keys name a prefix rather than a file; these are the
# suffixes the evaluation scripts append to it.
PREFIX_SUFFIXES: Tuple[str, ...] = (
    ".expressed.gtf", ".expressed_kept.gtf", ".excluded.gtf", ".reduced.gtf", ".reduced.db",
)

DEFAULT_CONFIG_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "configs")
DEFAULT_MANIFEST = os.path.join(DEFAULT_CONFIG_DIR, "data_manifest.tsv")

HASH_CHUNK = 8 * 1024 * 1024


def strip_column_spec(value: str) -> str:
    """Drop a trailing `:col` / `:col:col,col` selector from a path value."""
    if os.path.exists(value):
        return value
    head = value
    while ":" in head:
        head = head.rsplit(":", 1)[0]
        if os.path.exists(head):
            return head
    return value


def expand_prefix(path: str) -> List[str]:
    """Expand a reduced-db style prefix into the files that actually exist."""
    if os.path.isfile(path):
        return [path]
    found = [path + suffix for suffix in PREFIX_SUFFIXES if os.path.isfile(path + suffix)]
    return found


def collect_input_paths(config_files: Iterable[str]) -> Set[str]:
    """Return every existing input file referenced by the given configs."""
    paths: Set[str] = set()
    for config_file in config_files:
        try:
            data = yaml.safe_load(open(config_file))
        except Exception as e:  # noqa: BLE001 - a broken config should not abort the scan
            log.warning("Could not parse %s: %s", config_file, e)
            continue
        if not isinstance(data, dict):
            continue
        for key, value in data.items():
            if key in EXCLUDED_KEYS or not isinstance(value, str) or not value:
                continue
            if key not in SINGLE_PATH_KEYS and key not in MULTI_PATH_KEYS:
                continue
            tokens = value.split(" ") if key in MULTI_PATH_KEYS else [value]
            for token in tokens:
                if not token.startswith("/"):
                    continue
                paths.update(expand_prefix(strip_column_spec(token)))
    return paths


def sha256_of(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(HASH_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(manifest_file: str) -> Dict[str, Tuple[int, str]]:
    """Read a manifest into {path: (size, sha256)}."""
    entries: Dict[str, Tuple[int, str]] = {}
    if not os.path.exists(manifest_file):
        return entries
    for line in open(manifest_file):
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 3:
            continue
        entries[parts[0]] = (int(parts[1]), parts[2])
    return entries


def write_manifest(manifest_file: str, entries: Dict[str, Tuple[int, str]]) -> None:
    with open(manifest_file, "w") as f:
        f.write("# CI test input manifest - generated by data_manifest.py\n")
        f.write("# path\tsize_bytes\tsha256\n")
        for path in sorted(entries):
            size, digest = entries[path]
            f.write("%s\t%d\t%s\n" % (path, size, digest))


def config_files_in(config_dir: str) -> List[str]:
    return sorted(glob.glob(os.path.join(config_dir, "*.yaml")) +
                  glob.glob(os.path.join(config_dir, "*", "*.yaml")))


def do_generate(config_dir: str, manifest_file: str, reuse: bool) -> int:
    configs = config_files_in(config_dir)
    log.info("Scanning %d configs in %s", len(configs), config_dir)
    paths = collect_input_paths(configs)
    log.info("Found %d referenced input files", len(paths))

    known = load_manifest(manifest_file) if reuse else {}
    entries: Dict[str, Tuple[int, str]] = {}
    total = len(paths)
    hashed_bytes = 0
    for i, path in enumerate(sorted(paths), 1):
        size = os.path.getsize(path)
        cached = known.get(path)
        if cached and cached[0] == size:
            entries[path] = cached
            continue
        log.info("[%d/%d] hashing %s (%.1f MB)", i, total, path, size / 1048576)
        entries[path] = (size, sha256_of(path))
        hashed_bytes += size

    write_manifest(manifest_file, entries)
    log.info("Wrote %d entries to %s (%.1f GB hashed this run)",
             len(entries), manifest_file, hashed_bytes / 1073741824)
    return 0


def do_verify(config_dir: str, manifest_file: str, deep: bool,
              single_config: Optional[str]) -> int:
    manifest = load_manifest(manifest_file)
    if not manifest:
        log.warning("Manifest %s is empty or absent - nothing to verify", manifest_file)
        return 0

    if single_config:
        paths = collect_input_paths([single_config])
        scope = single_config
    else:
        paths = collect_input_paths(config_files_in(config_dir))
        scope = config_dir

    missing_file: List[str] = []
    size_mismatch: List[str] = []
    hash_mismatch: List[str] = []
    untracked: List[str] = []

    for path in sorted(paths):
        if path not in manifest:
            untracked.append(path)
            continue
        if not os.path.exists(path):
            missing_file.append(path)
            continue
        expected_size, expected_hash = manifest[path]
        if os.path.getsize(path) != expected_size:
            size_mismatch.append(path)
            continue
        if deep and sha256_of(path) != expected_hash:
            hash_mismatch.append(path)

    log.info("Verified %d inputs referenced by %s (%s)",
             len(paths), scope, "size + sha256" if deep else "size only")
    for label, items in (("MISSING", missing_file), ("SIZE CHANGED", size_mismatch),
                         ("HASH CHANGED", hash_mismatch)):
        for path in items:
            log.error("%s: %s", label, path)
    for path in untracked:
        log.warning("not in manifest (run `generate` to add): %s", path)

    if missing_file or size_mismatch or hash_mismatch:
        return IsoQuantExitCode.INPUT_FILE_NOT_FOUND
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=["generate", "verify"], help="what to do")
    parser.add_argument("--configs", type=str, default=DEFAULT_CONFIG_DIR,
                        help="directory holding the CI configs [%(default)s]")
    parser.add_argument("--manifest", type=str, default=DEFAULT_MANIFEST,
                        help="manifest file to write or check [%(default)s]")
    parser.add_argument("--config-file", type=str,
                        help="verify only the inputs of this single config")
    parser.add_argument("--deep", action="store_true",
                        help="verify by re-hashing, not just by size")
    parser.add_argument("--rehash-all", action="store_true",
                        help="generate: re-hash everything instead of reusing unchanged entries")
    return parser.parse_args()


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    args = parse_args()
    if args.mode == "generate":
        return do_generate(args.configs, args.manifest, reuse=not args.rehash_all)
    return do_verify(args.configs, args.manifest, args.deep, args.config_file)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(130)
