# Markdown summary of the result directories written by compare.jl:
#   julia benchmark/summary.jl benchmark/results/<run>
using TOML

med(xs) = (s = sort(Float64.(xs)); n = length(s); isodd(n) ? s[cld(n, 2)] : (s[n÷2] + s[n÷2+1]) / 2)

function tally(path)
    isfile(path) || return "–"
    rows = [split(l, '	') for l in readlines(path)[2:end]]
    st = last.(rows)
    e = [parse(Int, r[8]) for r in rows if r[end] in ("matched", "coordinates")]
    n(x) = count(==(x), st)
    "$(n("matched")) / $(n("coordinates")) / $(n("missing")) / $(n("extra"))" *
        (isempty(e) ? "" : " (≤$(maximum(e)) bp)")
end

root = isempty(ARGS) ? joinpath(@__DIR__, "results") : ARGS[1]
println("| run | RibosomalRNAFinder first / warm / forward (s) | nhmmer (s) | barrnap (s) | vs nhmmer: matched / endpoints >10 bp / missing / extra | vs barrnap |")
println("|---|---|---|---|---|---|")
for d in sort(readdir(root; join=true))
    f = joinpath(d, "timings.toml"); isfile(f) || continue
    t = TOML.parsefile(f)
    r(x) = string(round(x; digits=2))
    fw = haskey(t, "findrna_forward_warm_seconds") ? r(med(t["findrna_forward_warm_seconds"])) : "–"
    println("| ", basename(d), " | ", r(t["findrna_first_call_seconds"]), " / ", r(med(t["findrna_warm_seconds"])), " / ", fw,
            " | ", r(med(t["nhmmer_process_seconds"])), " | ", r(med(t["barrnap_process_seconds"])),
            " | ", tally(joinpath(d, "nhmmer-differences.tsv")), " | ", tally(joinpath(d, "barrnap-differences.tsv")), " |")
end
