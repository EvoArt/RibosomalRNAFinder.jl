# FindRNA.jl

Find ribosomal RNA genes in genomes and assemblies, as
[barrnap](https://github.com/tseemann/barrnap) does, from a plain Julia session:
no HMMER, Perl or WSL, so it works the same on Windows, macOS and Linux. It is
meant for the occasional handful of genomes in an interactive (e.g. VS Code)
session as much as for scripts, and fits alongside
[DenoiseDNA.jl](https://github.com/EvoArt/DenoiseDNA.jl) and
[AssignTaxonomy.jl](https://github.com/EvoArt/AssignTaxonomy.jl).

```julia
using FindRNA
hits = findrna("genome.fna")                     # 16S, 23S and 5S on both strands
```
```
RNAResult · 25 RNA hits
  source      genome.fna
  models      16S_rRNA × 8, 23S_rRNA × 8, 5S_rRNA × 9
  sequence    model         coordinates  strand    bits
  CP148752.1  16S_rRNA    111206:112737       +  1393.8
  CP148752.1  23S_rRNA    113330:116223       +  2556.1
  …
```

```julia
hits16 = findrna("genome.fna"; models=["16S_rRNA"])
seqs = extract_sequences(hits16)                 # Vector{LongDNA{4}}, 5'→3' along the RNA
write_fasta("16S.fasta", hits16)
write_gff("rrna.gff", hits)                      # barrnap-style GFF3

using DataFrames; DataFrame(hits)                # results are a Tables.jl table

using AssignTaxonomy                             # classify the 16S genes
assign_taxonomy(seqs, hit_ids(hits16), "rdp_train_set.fa")
```

Each `RNAHit` has `seqid`, `model`, `start`, `stop` (1-based, inclusive, on the
input sequence), `strand`, `bits`, `evalue` and `partial`.

## Options

Defaults follow barrnap 0.9:

| keyword     | default | meaning |
|-------------|---------|---------|
| `models`    | all     | profiles to search, e.g. `["16S_rRNA"]` |
| `evalue`    | `1e-6`  | report hits with E-value at most this |
| `reject`    | `0.25`  | drop hits shorter than this fraction of the model |
| `lencutoff` | `0.8`   | mark hits shorter than this fraction as `partial` |
| `profiles`  | bacterial 16S/23S/5S | any HMMER3 nucleotide profile file |
| `scoring`   | `:viterbi` | `:forward` rescores genes as nhmmer does (bits, bias and E-values within ~1 bit of nhmmer's; ~10% slower) |
| `threaded`  | `true`  | use all Julia threads (start Julia with `-t auto`) |
| `progress`  | interactive | progress bar; `FindRNA.setprogress!(false)` turns it off |

The bundled profiles are barrnap 0.9's bacterial models (from Rfam). Archaeal,
mitochondrial or eukaryotic HMMER3 profiles can be passed with `profiles`.

## Accuracy and speed

Searches use [ProfileHMMs.jl](https://github.com/EvoArt/ProfileHMMs.jl), a
native Julia implementation of the nhmmer pipeline: SIMD ungapped seeds, then
local Viterbi alignment around them. On five bacterial genomes (3.7–6.8 Mb,
GC 58–68%), FindRNA reports the same 69 genes as nhmmer 3.4 and barrnap 0.9
with the same profiles, with no extra or missing calls and endpoints within
6 bp. CI checks this on every push; see [benchmark/](benchmark/README.md).

A 6.8 Mb genome takes about 1 s with 4 threads (1.5 s on one), 4–6 times
faster than nhmmer with the same threads; the first call in a session
is barely slower, because the package is precompiled.

By default bit scores are Viterbi scores without nhmmer's composition-bias correction,
so they differ from nhmmer's by a few percent. With `scoring=:forward` FindRNA rescores
each gene exactly as nhmmer does, and bits and bias agree with nhmmer's to about a bit.
Modern barrnap (1.x) uses Infernal covariance models instead; its calls are
compared but not required to agree.

## Installation

```julia
using Pkg
Pkg.add("FindRNA")
```

ProfileHMMs is registered. Until FindRNA's first registration merges, install it
with `Pkg.add(url="https://github.com/EvoArt/FindRNA.jl")`. In a checkout of this
repository, run `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.

## Licence

Code: MIT. Rfam model data: CC0; see [data/README.md](data/README.md).
