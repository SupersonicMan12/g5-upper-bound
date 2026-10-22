# A certified upper bound $g_5 \le 4.15$ for complete pivoting

This repository contains the paper, node logs, certificate-generation code and an independent exact checker for the computer-assisted proof of

$$
g_5 \le 4.15,
$$

where $g_5$ is the largest growth factor of Gaussian elimination with complete pivoting
on real $5\times5$ matrices. More precisely, no such matrix, with any legal pivoting
sequence, has fifth pivot $p_5 \ge K$ for

- $K = 9/2 = 4.5$, and
- $K = 2175795/2^{19} = 4.1499996\ldots$ (hence $g_5 \le 4.15$).

The previous best upper bound was $g_5 \le 4.84$ (Chen, Edelman, Urschel, 2026); their
explicit matrix gives $g_5 \ge 4.1325\ldots$.

The certificates are too large for GitHub and are needed to run the checker. They are archived on Zenodo: [doi:10.5281/zenodo.23499516](https://doi.org/10.5281/zenodo.23499516).

## Quick verification of the result

Requirements: Python 3.9 or newer, `xz` and about 9 GB of free disk space (the downloaded parts can be deleted after joining). No other software is needed.

```sh
git clone https://github.com/SupersonicMan12/g5-upper-bound
cd g5-upper-bound/data
# download all files of https://doi.org/10.5281/zenodo.23499516 into data/:
#   certificates_K4p5.jsonl.xz, certificates_K4p15.jsonl.xz.part00 ... part49, SHA256SUMS
cat certificates_K4p15.jsonl.xz.part* > certificates_K4p15.jsonl.xz
xz -dk certificates_K4p5.jsonl.xz certificates_K4p15.jsonl.xz
shasum -a 256 -c SHA256SUMS # every line should print OK
cd ..
python3 checker.py 9/2 data/logs_K4p5 data/certificates_K4p5.jsonl
python3 checker.py 2175795/524288 data/logs_K4p15 data/certificates_K4p15.jsonl
```

On Zenodo, `certificates_K4p15.jsonl.xz` is stored as 50 consecutive byte chunks (`.part00` to `.part49`, about 20 MB each, made with `split -d -a 2 -n 50`); `cat` rejoins them in order, and `SHA256SUMS` confirms that the joined file is identical to the original. The `SHA256SUMS` on Zenodo is identical to the one in this repository.

On Linux, `sha256sum -c SHA256SUMS` does the same as `shasum -a 256 -c`. `xz` is included in most Linux distributions; on macOS it can be installed with Homebrew (`brew install xz`). On Windows, run the commands in Git Bash or WSL.

A successful run ends with "verified ... tightenings and ... leaves", and any failed check stops with an AssertionError. Measured on an Apple M4 with Python 3.11.9:

| $K$ | nodes | leaves | tightenings | certificates | time |
|---|---|---|---|---|---|
| $9/2$ | 7,543 | 3,772 | 99,127 | 102,899 | 51 s |
| $2175795/2^{19}$ | 225,465 | 112,733 | 3,366,236 | 3,478,969 | 32.5 min |

Full transcripts are in `results/`.

## Data

| file | size | SHA-256 |
|---|---|---|
| `certificates_K4p5.jsonl` | 154 MB | `265017a55f45e9889dd3594633623c261e91b11db256bbed2512451a2ff80ea7` |
| `certificates_K4p15.jsonl` | 6.18 GB | `2e3d57ec400c37bf80a32fb7f747e9bd27bb3bb6b5798b1339b2f7acfab91055` |
| `certificates_K4p5.jsonl.xz` | 26 MB | `c12086e1565a6775badda2d3197c15e5ae9580e9daf000001a5a4859698dda0c` |
| `certificates_K4p15.jsonl.xz` | 1.01 GB | `f796fb185a9a36769ed9360e150d7cf6c1ccf9271face0f7bb642dd9939edc44` |

`certificates_K4p5.jsonl.xz` is downloaded from Zenodo directly and `certificates_K4p15.jsonl.xz` is obtained by joining its 50 parts; the decompressed `.jsonl` files are exactly the files that `checker.py` verified. These hashes are also in `data/SHA256SUMS` and on the Zenodo record. The node logs (the search trees) are in `data/logs_K4p5/` and `data/logs_K4p15/`. All formats are described in `data_formats.md`.

## Overview of the code

`checker.py` is the independent exact verifier. `checklist.md` lists the checks it makes, as prescribed in the paper.

In `generators/search_v1/`, `lubb.jl` computes the relaxations and processes the nodes, and `lubb_search.jl` is the parallel search driver. Both trees were built with:

```sh
julia --startup-file=no --project=. generators/search_v1/lubb_search.jl --K=4.5 --workers=9 --budget=1000 --output=search.toml --logdir=logs
julia --startup-file=no --project=. generators/search_v1/lubb_search.jl --K=4.15 --workers=9 --budget=1000 --output=search.toml --logdir=logs
```

`generators/export_v3/` is a later revision of version 1, not an independent implementation. Its `lubb_export.jl` reads the node logs written by version 1, re-solves the LPs of each node, checks every certificate exactly and writes the certificate files used by `checker.py`. The exporter's SHA-256 is recorded in the header line of each certificate file.

```sh
julia --startup-file=no --project=. generators/export_v3/lubb_export.jl --logdir=data/logs_K4p5 --K=4.5 --output=certificates_K4p5.jsonl
julia --startup-file=no --project=. generators/export_v3/lubb_export.jl --logdir=data/logs_K4p15 --K=4.15 --output=certificates_K4p15.jsonl
```

All code in `generators/` is written in Julia and was run with Julia 1.12.7. All package versions, including HiGHS, are pinned in `Project.toml` and `Manifest.toml` (run `julia --project=. -e 'using Pkg; Pkg.instantiate()'` once). Some packages in `Project.toml` are left over from earlier experiments and are not used. The scripts start their worker processes with `--project` set to their own directory, so to rerun them, copy `Project.toml` and `Manifest.toml` into `generators/search_v1/` and `generators/export_v3/` first (the scripts originally lived next to these files). None of this is needed to run `checker.py`.

## Reports

`reports/` contains the reports of the original runs, unchanged. The search reports record the source hashes, options, timings and node counts. Please note:

- They refer to the logs by their original path `lubb_runs/logs_K*`, and the certificate headers by `4.5/logs_K4p5` and `4.15/logs_K4p15`. These are the same data as `data/logs_K*`.
- Their `statement` field says "g5 < K". The claim of this work is the one in the paper, $g_5 \le 4.15$.
- The replay reports are a second check of the trees by related AI-written code. That code (`lubb_replay.jl`, and for 4.15 an earlier version of `lubb.jl`) is not included, and the replays are not part of the proof.

## Repository layout

```
checker.py                  independent exact checker
checklist.md                specification of the checks
data_formats.md             node-log and certificate formats
data/                       node logs, SHA256SUMS (certificates: Zenodo)
generators/                 search (search_v1) and certificate export (export_v3)
results/                    checker transcripts with timings
reports/                    search and replay reports from the original runs
g5.pdf                      the paper
Project.toml, Manifest.toml pinned Julia environment
```

## Use of AI

This project was carried out with substantial help from AI coding agents (Claude Opus 5.5 on Devin, Cognition AI).

- **Method and code (AI):** the overall approach (the 24-parameter LU normal form,
  the linear relaxation of the bilinear terms, branch-and-bound with LP-duality
  certificates) and an initial draft of the proof were produced by AI. The AI also
  wrote the search and certificate-generation code in `generators/`, ran the searches
  (reports in `reports/`), and wrote `data_formats.md`.
- **Verification and paper (author):** I worked through the AI's proof draft until I understood it fully, then rewrote it as the paper, checking every step. I then independently wrote the exact checker `checker.py` and its specification `checklist.md`. AI was also used to explain background material, to review drafts and code, for proofreading and notation cleanup, and to check citations. All final claims and text were decided by me.

Importantly, the AI-written generators are not trusted: they only *propose* certificates. The final result rests on the mathematical argument in the paper, checked by the author, and on `checker.py`, written independently by the author, which re-verifies every certificate in exact rational arithmetic.

## License

Code is released under the MIT License (see `LICENSE`). The certificates and node logs are released under CC BY 4.0.
