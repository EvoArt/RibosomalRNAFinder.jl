#!/bin/bash
#SBATCH --job-name=findrna-bench
#SBATCH --partition=hmq
#SBATCH --exclude=athena-h03
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
# Speed and output comparison of RibosomalRNAFinder with HMMER 3.4 nhmmer, barrnap 0.9
# (same profiles) and modern barrnap (Infernal CMs; reported, not asserted).
# Submit from a checkout: sbatch benchmark/athena.sh
# Tools are built once under $TOOLS (default ~/findrna-tools).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
TOOLS="${TOOLS:-$HOME/findrna-tools}"
export PATH="$HOME/julia-1.12.7/bin:$PATH"
mkdir -p "$TOOLS"
if [ ! -x "$TOOLS/bin/nhmmer" ]; then
    curl -fL --retry 3 https://eddylab.org/software/hmmer/hmmer-3.4.tar.gz -o "$TOOLS/hmmer.tar.gz"
    echo "ca70d94fd0cf271bd7063423aabb116d42de533117343a9b27a65c17ff06fbf3  $TOOLS/hmmer.tar.gz" | sha256sum -c -
    tar -xzf "$TOOLS/hmmer.tar.gz" -C "$TOOLS"
    (cd "$TOOLS/hmmer-3.4" && ./configure --prefix="$TOOLS" && make -j4 && make install)
fi
[ -d "$TOOLS/barrnap-0.9" ] || git clone --depth 1 --branch 0.9 https://github.com/tseemann/barrnap.git "$TOOLS/barrnap-0.9"
if [ ! -x "$TOOLS/modern/bin/barrnap" ]; then
    "$HOME/miniforge3/bin/conda" create -y --override-channels -c conda-forge -c bioconda \
        -p "$TOOLS/modern" barrnap=1.10.6 hmmer=3.4
fi
if [ ! -x "$TOOLS/deps/bin/bedtools" ]; then     # barrnap 0.9 needs bedtools
    "$HOME/miniforge3/bin/conda" create -y --override-channels -c conda-forge -c bioconda -p "$TOOLS/deps" bedtools
fi
export NHMMER="$TOOLS/bin/nhmmer"
base_path="$PATH"
export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK" BENCH_THREADS="$SLURM_CPUS_PER_TASK" BENCH_REPEATS=3
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
benchmark/fetch_genomes.sh genomes
out="benchmark/results/athena-$SLURM_JOB_ID"
failed=0
# 4 threads for both barrnap versions, then single-threaded (0.9 only) for
# a one-core comparison with nhmmer --cpu 1.
for variant in 0.9 modern t1-0.9; do
    t=$SLURM_CPUS_PER_TASK; [ "$variant" = t1-0.9 ] && t=1
    export JULIA_NUM_THREADS=$t BENCH_THREADS=$t
    if [ "$variant" != modern ]; then
        export BARRNAP="$TOOLS/barrnap-0.9/bin/barrnap" BARRNAP_MODERN=false
        export PATH="$TOOLS/bin:$TOOLS/deps/bin:$base_path"
    else
        export BARRNAP="$TOOLS/modern/bin/barrnap" BARRNAP_MODERN=true
        export PATH="$TOOLS/modern/bin:$base_path"   # its own perl and Infernal
    fi
    FINDRNA_RESULTS="$out/$variant-synthetic" julia --startup-file=no --project=. benchmark/compare.jl || failed=1
    for genome in genomes/*.fasta; do
        FINDRNA_RESULTS="$out/$variant-$(basename "$genome" .fasta)" \
            julia --startup-file=no --project=. benchmark/compare.jl "$genome" || failed=1
    done
done
julia --startup-file=no --project=. benchmark/summary.jl "$out"
echo "PARITY_FAILURE=$failed"
exit "$failed"
