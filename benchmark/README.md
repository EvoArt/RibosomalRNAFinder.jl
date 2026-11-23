# Comparisons with HMMER and barrnap

`compare.jl` runs RibosomalRNAFinder, nhmmer 3.4 and barrnap on the same FASTA with the
same profiles (`data/bac.hmm`) and E-value cutoff (1e-5). It asserts that every
reference call has exactly one RibosomalRNAFinder call with the same sequence, model and
strand, endpoints within 10 bp, and that RibosomalRNAFinder makes no other calls.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. benchmark/compare.jl                  # synthetic fixture
benchmark/fetch_genomes.sh genomes                      # accessions in genomes.txt
julia --project=. benchmark/compare.jl genomes/CP148752.1.fasta
julia benchmark/summary.jl benchmark/results            # markdown table
```

Environment: `NHMMER`, `BARRNAP` (executables), `BENCH_THREADS`,
`BENCH_REPEATS`, `FINDRNA_RESULTS` (output directory). With
`BARRNAP_MODERN=true` and a barrnap 1.x executable, the barrnap calls are
reported (`barrnap-differences.tsv`) but asserted only on the synthetic
fixture: barrnap 1.x uses Infernal covariance models.

Each output directory keeps the raw nhmmer table, both GFFs, per-hit
differences, tool versions, input/model SHA-256 and timings. RibosomalRNAFinder's first
call (including any compilation) is reported separately from warm calls;
reference timings are whole-process wall times.

The synthetic fixture holds model-consensus sequences on both strands, two
copies in one contig, a repetitive negative, an all-N sequence, and a public
5S sequence (NCBI CP148752.1:116377-116486) inside a 5 kb window.

Latest Athena run: [athena-results.md](athena-results.md).

CI (`.github/workflows/CI.yml`) runs the fixture and the five genomes on every
push. `athena.sh` runs both barrnap versions on the University of Exeter Athena
cluster: `sbatch benchmark/athena.sh` from a checkout.

HMMER 3.4: <https://hmmer.org>; barrnap: <https://github.com/tseemann/barrnap>.
