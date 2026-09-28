# fARGene: Python 2.7 → Python 3 port — tracking notes

Original tool: [fannyhb/fargene](https://github.com/fannyhb/fargene) (Python 2.7).
Ported copy: `/Users/indadaz@qut.edu.au/Claude/fargenev2/` (this repo).
Untouched original: `/Users/indadaz@qut.edu.au/Claude/DeepARG_short-reads_test/fargene/` (still Python 2.7, unmodified).

## Why this was needed

- bioconda's published `fargene` recipe has `skip: True  # [py>=30]` — hardcoded to never build for Python 3.
- Separately, that recipe doesn't even pull from the real upstream repo: it sources from `thanhleviet/fargene` (an unrelated fork) pinned to that fork's stale `v0.1` tag, which ships only 11 of the 23 HMM models. It's missing all 9 aminoglycoside models (`aminoglycoside_model_a`–`i`) and all 3 macrolide models (`erm_typeA`, `erm_typeF`, `mph`).
- This port sources from real upstream (`fannyhb/fargene`, current `master`), so it has the complete 23-model set.

## Code changes: Python 2 → Python 3 (mechanical)

Applied across `fargene_analysis/*.py` and `fargene_model_creation/*.py`:

- `print` statements → `print()` calls
- `dict.has_key(x)` → `x in dict`
- `dict.iteritems()` → `dict.items()`
- `itertools.izip()` → `zip()`
- `raw_input()` → `input()`
- Implicit relative imports (`from utils import ...`) → explicit (`from .utils import ...`)
- `distutils.spawn.find_executable` → `shutil.which` (distutils was removed entirely in Python 3.12 — this import would hard-fail on any current interpreter)
- Fixed inconsistent tab/space indentation in `fargene_analysis.py` (Python 3 rejects mixed indentation that Python 2 tolerated)
- `setup.py`: added `python_requires='>=3.7'`

## Real bugs found and fixed during the port

1. **Multiprocessing pickling bug** (`fargene_model_creation/estimate_sensitivity.py`, `sort_one_hmmerfile()`): returned a `dict.values()` view, which is not picklable in Python 3 and would have crashed the moment `multiprocessing.Pool` tried to send it back to the parent process. Fixed by wrapping in `list()`.
2. **String-vs-numeric score comparison bug** (`fargene_analysis/utils.py`, `orf_classifier()`): compared HMM scores as raw strings instead of floats (`"9.5" < "10.2"` is `False` lexicographically — wrong). This bug existed in the original Python 2.7 code too; it just never manifested because the only reason it wasn't caught earlier is it only triggers when two candidate ORFs share the same sequence-region identifier, a case the standard test data never produced. Fixed by casting to `float` at parse time. Verified with a synthetic test exercising the exact failure case in both score orderings.
3. **`--no-orf-predict` was a no-op in metagenomic mode** (`fargene_analysis/fargene_analysis.py`, `parse_fastq_input()`): the ORF-prediction step after SPAdes assembly was never actually gated by `options.orf_predict`, only by whether a contig was retrieved. This was already true in the original Python 2.7 code; it was invisible before because ORFfinder (the only thing that step called) was never installed in any environment we tested, so the step always failed silently regardless of the flag. It only became a real, consequential bug once the prodigal fallback (below) made that step actually succeed by default. Fixed by adding the `options.orf_predict` gate.

## ORFfinder replaced with a prodigal fallback for metagenomic mode

The metagenomic pipeline's final ORF-calling step called NCBI's `ORFfinder` unconditionally, with no way to substitute another tool. `ORFfinder` isn't installable on this machine (NCBI only distributes a `linux-i64` binary now, no macOS build; this machine is Apple Silicon) — and it's also not distributed via any conda channel, so anyone without a Linux x86_64 box hits the same wall.

Investigation: prodigal's actual gene-finding was never the problem. Given the same assembled contig, prodigal correctly identified the exact same 816bp NDM-1 gene ORFfinder would have — but fARGene's own `parse_prodigal()` filter (`fargene_analysis/predict_orfs.py`) only accepted gene calls prodigal marked fully complete (`partial=00`), and short metagenomic contigs routinely get assembled right up to (or just past) the true gene boundary, so prodigal marks the call `partial` purely because it runs into the edge of the contig — not because the gene itself is incomplete. That's almost certainly why the original authors reached for ORFfinder specifically for this step instead of reusing prodigal (which they already use, by default, for genomic-mode input).

Fix: `parse_prodigal()` gained an `allow_partial` parameter. When enabled (used only for the metagenomic path, on by default there):
- gene calls truncated at only one end (`partial=10` or `partial=01`) are kept if prodigal's own confidence score (`conf=`) is at least **50** (out of 100) — chosen deliberately loose, since this step only proposes *candidate* sequences; the real filter is the resistance-gene HMM score threshold applied downstream, so favoring recall here costs nothing but a little wasted compute on candidates that get rejected later anyway
- calls truncated at **both** ends (`partial=11` — no confirmed start or stop codon at all) are always rejected as too speculative
- fully complete calls (`partial=00`) are accepted exactly as before

`--orf-finder` still exists and now works consistently in both genomic and metagenomic mode: pass it to force the old ORFfinder-based behavior on a machine that actually has it installed (e.g. Linux). Default (no flag) is prodigal in both modes.

Verified on the tutorial's paired-end dataset: the metagenomic pipeline now runs fully to completion without ORFfinder and recovers an 819bp ORF that contains the exact same 813bp NDM-1 coding sequence found via genomic mode, byte-for-byte, with 6 extra bases at the 5' end (an alternate `TTG` start codon prodigal chose when working from a short, context-free contig fragment rather than a complete genome — a normal, expected variation in de novo start-codon calling, not a correctness problem).

This is a deliberate, disclosed departure from "same output as Python 2.7" for the metagenomic ORF step specifically — the original code path (ORFfinder-only) remains ported and available via `--orf-finder`, just unverifiable on this machine.

## Optimizations (all re-verified to produce byte-identical output to the pre-optimization port)

1. `orf_classifier()` (`fargene_analysis/utils.py`): was O(n²) — rescanned every existing hit on every new hit to find a same-region match. Replaced with an O(n) index dict.
2. `create_subsets()` (`fargene_model_creation/estimate_sensitivity.py`): was re-reading the entire reference FASTA from disk once per reference sequence (O(n²) I/O) to build each leave-one-out subset. Fixed to read the file once into memory.
3. Leave-one-out model building (`fargene_model_creation/estimate_sensitivity.py`): was fully serial (align → build → search, one model at a time). Parallelized with `multiprocessing.Pool`, mirroring the pooling pattern already used elsewhere in the same file. Measured: **161s → 102s** (~37% faster) on the tutorial dataset.
4. Removed shell pipeline overhead in three hot-path functions (`fargene_analysis/utils.py`): `translate_and_search()`, `classifier()`, and `orf_classifier()`'s hit-file generation previously shelled out through `cat | transeq | hmmsearch` and `grep | awk | awk` for every input file. Replaced with direct subprocess chaining / native Python filtering. Measured: genomic golden-path run dropped from ~15s to **2.3s**.
5. Minor: replaced a `wc -l` subprocess spawn (`ResultsSummary.count_hits()`) with a pure-Python line count; removed a stray leftover debug `print(logger.handlers)`; removed dead/unreachable debug code (`calculate_performance.py`'s `__main__` block referenced an undefined `Estimator` class and hardcoded paths from the original author's machine — never worked, never executed in the real pipeline).

## Software / environment changes

- **New conda environment**: `fargene_env_py3` (`/Users/indadaz@qut.edu.au/miniconda3/envs/fargene_env_py3`) — Python 3.10, numpy, matplotlib, hmmer 3.4, emboss 6.6.0, seqtk, spades, trim-galore, prodigal, clustalo. The ported package is `pip install -e`'d into it from `fargenev2/`.
- **Existing `fargene_env` (Python 2.7) — one change made**: installed `clustalo` into it (via `conda install -n fargene_env -c bioconda -c conda-forge clustalo`). It was missing entirely, which meant `fargene_model_creation` had never actually been runnable in that environment before this session, independent of anything else. This was needed to get a true Python 2.7 reference run for side-by-side comparison. No other changes were made to that environment or to the original `fargene/` code.
- **`environment.yml`** in `fargenev2/` updated to describe the Python 3 env (bumped `python=2.7` → `python=3.10`, added `clustalo` for `fargene_model_creation`, kept the existing hmmer/emboss version-pin annotations from earlier debugging).
- **ORFfinder**: not installed anywhere, and could not be made to work on this machine — NCBI now only distributes a `linux-i64` binary (no macOS build published anymore), and this machine is Apple Silicon (`osx-arm64`), so it wouldn't run natively even if fetched. The ORFfinder code path itself (`--orf-finder` flag) is ported and compiles cleanly but is **not runtime-verified**. However, the metagenomic pipeline no longer depends on it by default — see "ORFfinder replaced with a prodigal fallback" above.
- **pixi**: installed (`curl -fsSL https://pixi.sh/install.sh | bash`) to build and validate `pixi.toml`/`pixi.lock` in this repo, which give a one-command, fully reproducible install (`pixi install`) across Linux, Intel macOS, and Apple Silicon macOS — resolving both the conda-level prerequisites and the Python package itself. `pixi.lock` was actually solved (not hand-written) via `pixi install` and verified end-to-end with `pixi run test-tutorial` (same tutorial golden path as elsewhere in this document), which reproduced the expected 1 predicted gene / 1 predicted ORF result.

## Verification performed

All comparisons below are byte-for-byte diffs against genuine Python 2.7 runs (not just against the tutorial's documented reference numbers), using the tool's own tutorial dataset:

| Path | Result |
|---|---|
| Genomic mode (`--hmm-model class_b_1_2`, single FASTA) | byte-for-byte identical output files, before and after optimization |
| Metagenomic mode (paired-end FASTQ, multiprocessing pool, `Transformer`, quality trimming, SPAdes assembly) | identical through every stage; final assembly stage validated against the tutorial's documented reference value (the actual old-env SPAdes 3.15.2 binary crashes on this Mac for unrelated reasons — a pre-existing, pre-session environment issue, not a porting artifact) |
| `fargene_model_creation` (clustalo/hmmbuild/hmmpress/leave-one-out cross-validation, incl. the multiprocessing pickling fix) | deterministic (full-length) results byte-for-byte identical; fragment-based results match within expected random-sampling noise (no fixed seed in the original code, so exact reproduction isn't expected even between two Python 2.7 runs) |
| `--protein` flag (real peptide input) | byte-for-byte identical |
| `--store-peptides` + `--rerun` flags | byte-for-byte identical |
| Metagenomic ORF prediction (prodigal fallback, default) | full pipeline run-to-completion; recovered ORF contains the exact same 813bp gene sequence as genomic mode, byte-for-byte (see above for the 6bp 5'-end difference and why) |
| `--orf-finder` flag (ORFfinder path) | ported, compiles, exercised enough to confirm it still correctly requires/requests ORFfinder and fails gracefully when absent — **not runtime-tested with a working ORFfinder binary** (platform limitation, see above) |

## Upstream repo status (`fannyhb/fargene`)

Checked directly against the GitHub repo and API:

- Two existing tags: `v0.1` and `v1.0`. **Neither has the full 23-model set** — both only ship the same 11 models (beta-lactamases, qnr, tet genes). Only current `master` (commit `41e3386`) has the aminoglycoside and macrolide models added.
- Repo last pushed **2022-07-14** — dormant for 4+ years as of this writing.
- Open, unanswered issues asking for exactly what this port addresses: [#18 "Request for new release"](https://github.com/fannyhb/fargene/issues/18), [#16 "Missing pre-defined models"](https://github.com/fannyhb/fargene/issues/16), [#19 "conda installation seems it works but not"](https://github.com/fannyhb/fargene/issues/19) — all opened by other users independently hitting the same problems, none acted on.

Conclusion: waiting on upstream to cut a new release isn't realistic based on the evidence above. The practical path is to publish this Python 3 port (with the model-set fix, the two bug fixes, and the ORFfinder→prodigal fallback) as its own citable release — either as a PR back to `fannyhb/fargene` (worth attempting once, given the repo isn't archived) or, if unanswered, as a standalone tagged release under a new home — and point the bioconda recipe at *that*.

## What "fix the bioconda recipe" actually requires

1. **A citable, tagged source to point at.** Bioconda recipes pin to an immutable release archive (`https://github.com/<org>/<repo>/archive/v{{ version }}.tar.gz` + a `sha256`), not a moving branch. This port doesn't have one yet — the first concrete step is publishing this `fargenev2/` code (once you're happy with it) as a tagged GitHub release somewhere citable.
2. **Fork `bioconda/bioconda-recipes`** on GitHub, clone it, branch off `master`.
3. **Edit `recipes/fargene/meta.yaml`**:
   - Remove `build: skip: True  # [py>=30]`
   - Change `source: url:` from `thanhleviet/fargene`'s stale `v0.1` fork tag to the new tagged release from step 1, and update the `sha256` to match the new tarball (`sha256sum` on the downloaded archive)
   - Bump `package: version:` and/or `build: number:` to reflect the new source
   - Review `requirements: run:` — the current list includes `biopython` and `pyyaml`, neither of which this codebase actually imports anywhere (confirmed by grepping all `import`/`from` statements) — worth dropping as dead weight, though not required for correctness
   - `ORFfinder` is currently fetched only for `# [linux]` in the `source:` section (no macOS build exists from NCBI) — worth flagging in the recipe or docs that `--orf-finder` is Linux-only, now that macOS/other platforms have the prodigal fallback as the default
4. **Test locally** before opening a PR — `bioconda-utils lint` and `bioconda-utils build` against the changed recipe (bioconda's own CI runs these; catching issues locally first avoids review churn)
5. **Open the PR** against `bioconda/bioconda-recipes`; it triggers bioconda's CI (linux-64 + osx-64 builds); address any review feedback from bioconda maintainers
6. Once merged, `conda install -c bioconda fargene` gets the fixed package for everyone

Not yet done: step 1 (no tagged release of this port exists yet) blocks everything after it. This is naturally sequenced *after* the paper, per your original plan.
