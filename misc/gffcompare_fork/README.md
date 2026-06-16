# gffcompare fork: `--terminal-delta`

`terminal_delta.patch` adds a `--terminal-delta <d>` option to gffcompare that
enforces a terminal-end tolerance `d` (bp) at the **transcript level only**,
independently of `-e` (which stays tied to exon-level Sn/Sp) and without the
other `--strict-match` side effects. It makes Transcript-level Sn/Sp
end-sensitive while Intron-chain level stays end-agnostic — needed so the
transcript-discovery benchmark can see polyA/TSS terminal-coordinate accuracy
(default gffcompare ignores terminal ends for multi-exon transcripts). See
`.claude/GFFCOMPARE.md` for the full rationale and the matching internals.

## Reproduce the binary

```bash
git clone --branch v0.12.6 https://github.com/gpertea/gffcompare.git gffcompare_term_delta
cd gffcompare_term_delta
git clone https://github.com/gpertea/gclib.git gclib   # Makefile clones gclib master
git apply /path/to/IsoQuant2/misc/gffcompare_fork/terminal_delta.patch
make GCLIB=gclib release
# install so the assessment pipeline (PATH) uses it:
cp gffcompare "$CI_BIN/gffcompare"     # back up the stock binary first
```

The patch touches only `gffcompare.cpp`, `gtf_tracking.cpp`, `gtf_tracking.h`
(a `terminalStrict` flag, the option parsing, gating the two `~`->`=`
promotions, and help text). The default behavior with no `--terminal-delta`
is byte-identical to stock v0.12.6.

## Usage

```
gffcompare --terminal-delta=10 -r <ref.gtf> -o <prefix> <query.gtf>
```

`misc/reduced_db_gffcompare.py` runs each split at deltas
`[None, 50, 10]` (default + the two end-sensitive levels) and
`isoquant_tests/github/run_pipeline.py` checks the three baseline blocks
`transcripts` / `transcripts_td50` / `transcripts_td10` in the CI YAML configs.
