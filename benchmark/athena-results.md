# Athena results (job 4052411, 2026-11-15)

Skylake-AVX512 node; Julia 1.12.7; HMMER 3.4; barrnap 0.9 and 1.10.6. Rows `0.9-*`
and `modern-*` use 4 threads for every tool, `t1-*` one thread. All tools
searched with E-value 1e-5; nhmmer and barrnap 0.9 used the same profiles as
RibosomalRNAFinder. Times are medians of 3 runs (whole-process wall time for nhmmer and
barrnap). "forward" is RibosomalRNAFinder with `scoring=:forward` (nhmmer-style scores);
its bits and bias were within 2 bits of nhmmer's for every gene
(`nhmmer-scores.tsv` in each result directory).

Rows marked `modern` compare against barrnap 1.10.6 (Infernal covariance
models); its 23S calls are absent for three genomes and its 5S boundaries
differ by up to 12 bp, so these rows are reported, not asserted.

| run | RibosomalRNAFinder first / warm / forward (s) | nhmmer (s) | barrnap (s) | vs nhmmer: matched / endpoints >10 bp / missing / extra | vs barrnap |
|---|---|---|---|---|---|
| 0.9-CP148752.1 | 1.16 / 0.94 / 1.01 | 4.53 | 3.99 | 25 / 0 / 0 / 0 (≤2 bp) | 25 / 0 / 0 / 0 (≤2 bp) |
| 0.9-CP148753.1 | 0.94 / 0.77 / 0.8 | 3.62 | 3.4 | 13 / 0 / 0 / 0 (≤3 bp) | 13 / 0 / 0 / 0 (≤3 bp) |
| 0.9-JBBHKO010000001.1 | 0.89 / 0.67 / 0.7 | 2.81 | 2.83 | 9 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 0 (≤3 bp) |
| 0.9-NZ_JBBHKP010000001.1 | 0.74 / 0.49 / 0.51 | 2.12 | 2.07 | 10 / 0 / 0 / 0 (≤6 bp) | 10 / 0 / 0 / 0 (≤5 bp) |
| 0.9-NZ_JBBHKQ010000001.1 | 0.77 / 0.53 / 0.57 | 2.55 | 2.24 | 12 / 0 / 0 / 0 (≤4 bp) | 12 / 0 / 0 / 0 (≤4 bp) |
| 0.9-synthetic | 0.36 / 0.11 / 0.12 | 0.72 | 0.85 | 9 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 0 (≤3 bp) |
| modern-CP148752.1 | 1.14 / 0.92 / 1.01 | 4.53 | 16.33 | 25 / 0 / 0 / 0 (≤2 bp) | 25 / 0 / 0 / 0 (≤4 bp) |
| modern-CP148753.1 | 0.96 / 0.75 / 0.84 | 3.59 | 11.13 | 13 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 4 (≤4 bp) |
| modern-JBBHKO010000001.1 | 0.87 / 0.66 / 0.7 | 2.88 | 10.32 | 9 / 0 / 0 / 0 (≤3 bp) | 6 / 0 / 0 / 3 (≤7 bp) |
| modern-NZ_JBBHKP010000001.1 | 0.73 / 0.5 / 0.53 | 2.11 | 8.94 | 10 / 0 / 0 / 0 (≤6 bp) | 6 / 4 / 0 / 0 (≤12 bp) |
| modern-NZ_JBBHKQ010000001.1 | 0.75 / 0.51 / 0.56 | 2.52 | 10.13 | 12 / 0 / 0 / 0 (≤4 bp) | 8 / 0 / 0 / 4 (≤4 bp) |
| modern-synthetic | 0.39 / 0.09 / 0.18 | 0.73 | 27.15 | 9 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 0 (≤4 bp) |
| t1-0.9-CP148752.1 | 1.92 / 1.53 / 1.77 | 10.3 | 8.83 | 25 / 0 / 0 / 0 (≤2 bp) | 25 / 0 / 0 / 0 (≤2 bp) |
| t1-0.9-CP148753.1 | 1.61 / 1.38 / 1.53 | 7.93 | 7.22 | 13 / 0 / 0 / 0 (≤3 bp) | 13 / 0 / 0 / 0 (≤3 bp) |
| t1-0.9-JBBHKO010000001.1 | 1.4 / 1.25 / 1.28 | 6.69 | 6.11 | 9 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 0 (≤3 bp) |
| t1-0.9-NZ_JBBHKP010000001.1 | 1.08 / 0.74 / 0.84 | 4.57 | 4.43 | 10 / 0 / 0 / 0 (≤6 bp) | 10 / 0 / 0 / 0 (≤5 bp) |
| t1-0.9-NZ_JBBHKQ010000001.1 | 1.11 / 0.93 / 1.05 | 5.25 | 4.72 | 12 / 0 / 0 / 0 (≤4 bp) | 12 / 0 / 0 / 0 (≤4 bp) |
| t1-0.9-synthetic | 0.24 / 0.11 / 0.19 | 0.72 | 0.87 | 9 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 0 (≤3 bp) |
