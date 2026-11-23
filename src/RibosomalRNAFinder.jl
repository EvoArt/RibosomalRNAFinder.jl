"""
    RibosomalRNAFinder

Find ribosomal RNA genes in FASTA files, in the manner of barrnap, from a
native Julia session on Windows, macOS or Linux. See [`findrna`](@ref).
"""
module RibosomalRNAFinder

using BioSequences, FASTX, Tables
using Printf: @sprintf
using ProfileHMMs
using ProfileHMMs: setprogress!
using PrecompileTools: @setup_workload, @compile_workload
import ProfileHMMs: _title, _kv, _table, _more, _nshow, _g, _fmt

export findrna, RNAHit, RNAResult, extract_sequences, hit_ids, write_fasta, write_gff
export default_profiles, setprogress!

"""
    RNAHit

One predicted RNA. `start:stop` are 1-based inclusive coordinates on the input
sequence `seqid`; `strand` is `'+'` or `'-'`. `partial` marks a hit shorter
than `lencutoff` of the model length. `bias` is the composition-bias correction
already subtracted from `bits` (computed only with `scoring=:forward`; else 0).
"""
struct RNAHit
    seqid::String
    model::String
    start::Int
    stop::Int
    strand::Char
    bits::Float64
    bias::Float64
    evalue::Float64
    partial::Bool
end
Base.length(h::RNAHit) = h.stop - h.start + 1
Base.show(io::IO, h::RNAHit) = print(io, "RNAHit(", h.model, ", ", h.seqid, ":", h.start, "-",
                                     h.stop, "(", h.strand, ")", h.partial ? ", partial" : "", ")")

"Hits with vector indexing and Tables.jl row access, plus the FASTA they came from."
struct RNAResult <: AbstractVector{RNAHit}
    hits::Vector{RNAHit}
    source::String
end
Base.size(r::RNAResult) = size(r.hits)
Base.getindex(r::RNAResult, i::Int) = r.hits[i]
Base.getindex(r::RNAResult, i::AbstractVector) = RNAResult(r.hits[i], r.source)
Base.IndexStyle(::Type{RNAResult}) = IndexLinear()
Base.filter(f, r::RNAResult) = RNAResult(filter(f, r.hits), r.source)
Tables.istable(::Type{RNAResult}) = true
Tables.rowaccess(::Type{RNAResult}) = true
Tables.rows(r::RNAResult) = r.hits
Tables.schema(::RNAResult) = Tables.Schema(fieldnames(RNAHit), fieldtypes(RNAHit))

"Path to the bundled bacterial 16S, 23S and 5S rRNA profiles (Rfam, as used by barrnap 0.9)."
default_profiles() = joinpath(@__DIR__, "..", "data", "bac.hmm")

"""
    findrna(fasta; models=nothing, evalue=1e-6, reject=0.25, lencutoff=0.8,
            profiles=default_profiles(), threaded=true, progress)

Find rRNA genes on both strands of every sequence in `fasta`.

- `models`: names of the profiles to search, e.g. `["16S_rRNA"]` (default: all).
- `evalue`: report hits with E-value at most this (barrnap's `--evalue`).
- `reject`: drop hits shorter than this fraction of the model (`--reject`).
- `lencutoff`: mark hits shorter than this fraction as `partial` (`--lencutoff`).
- `profiles`: a HMMER3 nucleotide profile file, for other RNA families.
- `scoring`: `:viterbi` (default, fastest) scores each gene by its best
  alignment. `:forward` rescores it as nhmmer (and so barrnap 0.9) does:
  Forward score with composition-bias correction and Forward-calibrated
  E-values, giving bits and E-values within about a bit of nhmmer's.

Further keywords (`seed_bits`, `db_size`) go to `ProfileHMMs.scan_fasta`.
Returns an [`RNAResult`](@ref): a vector of [`RNAHit`](@ref) and a Tables.jl table.
"""
function findrna(fasta::AbstractString; profiles=default_profiles(), models=nothing,
                 evalue::Real=1e-6, reject::Real=0.25, lencutoff::Real=0.8, kwargs...)
    0 <= reject <= 1 || throw(ArgumentError("reject must be in [0, 1]"))
    0 <= lencutoff <= 1 || throw(ArgumentError("lencutoff must be in [0, 1]"))
    isfile(fasta) || throw(ArgumentError("no such file: $fasta"))
    lengths = Dict(p.name => p.length for p in parse_hmm_file(profiles))
    hits = RNAHit[]
    for h in nhmmscan(profiles, fasta; models=models, evalue_threshold=evalue, kwargs...)
        M = lengths[h.query_name]
        len = h.target_end - h.target_start + 1
        len < floor(Int, reject * M) && continue
        push!(hits, RNAHit(h.target_name, h.query_name, h.target_start, h.target_end,
                           h.strand, h.bit_score, h.bias, h.evalue, len < floor(Int, lencutoff * M)))
    end
    RNAResult(hits, abspath(fasta))
end

"""
    extract_sequences(r::RNAResult; fasta=r.source) -> Vector{LongDNA{4}}

The sequence of each hit, in result order, oriented 5'→3' along the RNA
(reverse-complemented for `-` strand hits).
"""
function extract_sequences(r::RNAResult; fasta=r.source)
    wanted = Set(h.seqid for h in r)
    seqs = Dict{String,LongDNA{4}}()
    isempty(wanted) && return LongDNA{4}[]
    open(FASTA.Reader, fasta) do reader
        for record in reader
            id = String(identifier(record))
            id in wanted || continue
            haskey(seqs, id) && throw(ArgumentError("duplicate FASTA identifier: $id"))
            seqs[id] = LongDNA{4}(sequence(record))
        end
    end
    map(r.hits) do h
        haskey(seqs, h.seqid) || throw(ArgumentError("sequence $(h.seqid) not found in $fasta"))
        s = seqs[h.seqid][h.start:h.stop]
        h.strand == '-' ? reverse_complement(s) : s
    end
end

"""
    hit_ids(r::RNAResult) -> Vector{String}

An identifier per hit, `seqid:start-stop(strand)`, as written by
[`write_fasta`](@ref). With [`extract_sequences`](@ref) this gives the
`(sequences, ids)` pair taken by `AssignTaxonomy.assign_taxonomy`.
"""
hit_ids(r::RNAResult) = ["$(h.seqid):$(h.start)-$(h.stop)($(h.strand))" for h in r]

"Write the hit sequences as FASTA; headers are `hit_ids` followed by the model name."
function write_fasta(io::IO, r::RNAResult; kwargs...)
    for (id, h, s) in zip(hit_ids(r), r, extract_sequences(r; kwargs...))
        println(io, ">", id, " ", h.model, h.partial ? " (partial)" : "")
        println(io, s)
    end
    nothing
end

# GFF3 escaping as %XX: seqids allow only a small set of characters; attribute
# values must escape only separators, '%' and control characters.
_pct(b) = "%" * uppercase(string(b; base=16, pad=2))
const _SEQID_SAFE = Set(codeunits("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:^*\$@!+_?-|"))
_gff_seqid(s) = join(b in _SEQID_SAFE ? string(Char(b)) : _pct(b) for b in codeunits(s))
_gff_attr(s) = join(b in codeunits(";=&,%") || b < 0x20 || b == 0x7f ? _pct(b) : String([b])
                    for b in codeunits(s))

_product(model) = (m = match(r"^(\d+S)_rRNA$", model)) === nothing ? model : "$(m[1]) ribosomal RNA"

"Write GFF3 in barrnap's layout: column 6 is the E-value, `product` notes partial hits."
function write_gff(io::IO, r::RNAResult)
    println(io, "##gff-version 3")
    for h in r
        feature = endswith(h.model, "rRNA") ? "rRNA" : "ncRNA"
        attrs = "Name=$(_gff_attr(h.model));product=$(_gff_attr(_product(h.model)))" *
                (h.partial ? " (partial)" : "")
        println(io, join((_gff_seqid(h.seqid), "RibosomalRNAFinder", feature, h.start, h.stop,
                          h.evalue, h.strand, '.', attrs), '\t'))
    end
    nothing
end

for f in (:write_fasta, :write_gff)
    @eval function $f(path::AbstractString, r::RNAResult; kwargs...)
        open(path, "w") do io
            $f(io, r; kwargs...)
        end
        path
    end
end

Base.show(io::IO, r::RNAResult) = print(io, "RNAResult(", length(r), " hits)")
function Base.show(io::IO, ::MIME"text/plain", r::RNAResult)
    _title(io, "RNAResult", "$(_fmt(length(r))) RNA hits")
    _kv(io, "source", basename(r.source))
    isempty(r) && return
    models = unique(h.model for h in r)
    _kv(io, "models", join(("$m × $(count(h -> h.model == m, r))" for m in models), ", "))
    n = _nshow(io, length(r), 8)
    rows = [[first(h.seqid, 24), h.model, "$(h.start):$(h.stop)", string(h.strand),
             @sprintf("%.1f", h.bits), h.partial ? "partial" : ""] for h in r.hits[1:n]]
    _table(io, ["sequence", "model", "coordinates", "strand", "bits", ""], rows;
           align=[:l, :l, :r, :r, :r, :l])
    _more(io, length(r), n, "hits")
end

# The first `findrna` of a session should not wait for code generation.
@setup_workload begin
    p5 = last(parse_hmm_file(default_profiles()))
    cons = join("ACGT"[argmin(p5.match_emissions[:, k])] for k in 1:p5.length)
    @compile_workload begin
        path, io = mktemp()
        println(io, ">s\n", "N"^50, cons, "N"^50); close(io)
        r = findrna(path; progress=false)
        sprint(show, MIME"text/plain"(), r)
        sprint(write_gff, r); sprint(write_fasta, r)
        rm(path; force=true)
    end
end

end
