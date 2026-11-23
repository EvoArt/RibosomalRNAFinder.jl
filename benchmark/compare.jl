# Required references are never silently skipped. Output files retain raw evidence.
using RibosomalRNAFinder, ProfileHMMs, Test, SHA, TOML
include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

function parse_nhmmer(path)
    rows = NamedTuple[]
    for line in eachline(path)
        (startswith(line,"#") || isempty(strip(line))) && continue
        p = split(line)
        a,b = parse.(Int,p[7:8])
        push!(rows,(target=String(p[1]), model=String(p[3]), strand=only(p[12]),
                    start=min(a,b), stop=max(a,b), bits=parse(Float64,p[14]), bias=parse(Float64,p[15])))
    end
    rows
end

function parse_gff(path)
    rows = NamedTuple[]
    for line in eachline(path)
        (startswith(line,"#") || isempty(strip(line))) && continue
        p = split(line,'\t')
        p[3] == "rRNA" || continue
        m = match(r"(?:^|;)Name=([^;]+)",p[9])
        m === nothing && error("GFF record lacks model Name: $line")
        push!(rows,(target=String(p[1]),model=String(m[1]),strand=only(p[7]),
                    start=parse(Int,p[4]),stop=parse(Int,p[5]),bits=NaN))
    end
    rows
end

function compare_hits(got, ref; coordinate_tolerance=10)
    # One-to-one matching: an extra call fails just as a missing call does.
    unmatched = collect(eachindex(got))
    deltas = Int[]
    for r in ref
        candidates = filter(unmatched) do i
            h = got[i]
            h.seqid == r.target && h.model == r.model && h.strand == r.strand
        end
        @test !isempty(candidates)
        isempty(candidates) && continue
        errors = [max(abs(got[i].start-r.start),abs(got[i].stop-r.stop)) for i in candidates]
        k = argmin(errors); i = candidates[k]
        @test errors[k] <= coordinate_tolerance
        push!(deltas,errors[k]); filter!(!=(i),unmatched)
    end
    @test isempty(unmatched)
    deltas
end

function write_differences(path, got, refs)
    remaining = collect(eachindex(got))
    open(path,"w") do io
        println(io,"target\tmodel\tstrand\treference_start\treference_end\tjulia_start\tjulia_end\tendpoint_error\tjulia_bits\treference_bits\tstatus")
        for r in refs
            candidates = filter(i -> (got[i].seqid,got[i].model,got[i].strand) ==
                                     (r.target,r.model,r.strand),remaining)
            if isempty(candidates)
                println(io,join((r.target,r.model,r.strand,r.start,r.stop,"","","","",r.bits,"missing"),'\t'))
                continue
            end
            errs = [max(abs(got[i].start-r.start),abs(got[i].stop-r.stop)) for i in candidates]
            k = argmin(errs); i = candidates[k]; h = got[i]
            filter!(!=(i),remaining)
            println(io,join((r.target,r.model,r.strand,r.start,r.stop,h.start,h.stop,
                            errs[k],h.bits,r.bits,errs[k]<=10 ? "matched" : "coordinates"),'\t'))
        end
        for i in remaining
            h = got[i]
            println(io,join((h.seqid,h.model,h.strand,"","",h.start,h.stop,
                            "",h.bits,"","extra"),'\t'))
        end
    end
end

# Forward-mode bits and bias against nhmmer's, per matched hit: they use the
# same scoring, so they should agree to within about a bit.
function matched_scores(got, refs)
    map(refs) do r
        i = findfirst(h -> (h.seqid, h.model, h.strand) == (r.target, r.model, r.strand) &&
                           abs(h.start - r.start) <= 10 && abs(h.stop - r.stop) <= 10, got)
        (r, i === nothing ? nothing : got[i])
    end
end
function compare_scores(got, refs; tolerance=2.0)
    for (r, h) in matched_scores(got, refs)
        @test h !== nothing
        h === nothing && continue
        @test abs(h.bits - r.bits) <= tolerance
        @test abs(h.bias - r.bias) <= tolerance
    end
end
function write_scores(path, got, refs)
    open(path, "w") do io
        println(io, "target	model	strand	start	nhmmer_bits	findrna_bits	nhmmer_bias	findrna_bias")
        for (r, h) in matched_scores(got, refs)
            println(io, join((r.target, r.model, r.strand, r.start, r.bits, h === nothing ? "" : round(h.bits, digits=1),
                              r.bias, h === nothing ? "" : round(h.bias, digits=1)), '	'))
        end
    end
end

function median_time(xs)
    s = sort(xs); n = length(s)
    isodd(n) ? s[cld(n,2)] : (s[n÷2]+s[n÷2+1])/2
end

function benchmark_main()
    out = abspath(get(ENV,"FINDRNA_RESULTS",joinpath(@__DIR__,"results")))
    mkpath(out)
    fasta = isempty(ARGS) ? joinpath(out,"fixture.fna") : abspath(ARGS[1])
    isempty(ARGS) && make_fixture(fasta)
    hmmer = get(ENV,"NHMMER","nhmmer")
    barrnap = get(ENV,"BARRNAP","barrnap")
    threads = parse(Int,get(ENV,"BENCH_THREADS","1"))
    repeats = parse(Int,get(ENV,"BENCH_REPEATS","3"))
    repeats > 0 || error("BENCH_REPEATS must be positive")
    hv = read(`$hmmer -h`,String)
    occursin("HMMER 3.4",hv) || error("HMMER 3.4 required; got: $hv")
    version_path = joinpath(out,"barrnap-version.txt")
    open(version_path,"w") do io
        run(pipeline(`$barrnap --version`; stdout=io,stderr=io))
    end
    bv = read(version_path,String)
    write(joinpath(out,"versions.txt"),"Julia $VERSION\n$hv\n$bv")
    tbl = joinpath(out,"nhmmer.tbl")
    gff = joinpath(out,"barrnap.gff")
    hmm = default_profiles()
    # Equal input and profiles for HMMER/0.9. Modern Barrnap uses its own CMs.
    nhcmd = `$hmmer --cpu $threads --noali -E 1e-5 --tblout $tbl $hmm $fasta`
    bacmd = `$barrnap --threads $threads --kingdom bac --evalue 1e-5 $fasta`
    function external(cmd, output, err)
        run(pipeline(cmd; stdout=output, stderr=err))
    end
    opts = (progress=false,threaded=threads>1,evalue=1e-5)
    cold = @elapsed hits = findrna(fasta;opts...)
    println("RibosomalRNAFinder first call: ",cold," seconds, ",length(hits)," hits")
    flush(stdout)
    write_gff(joinpath(out,"findrna.gff"),hits)
    jt = [@elapsed(findrna(fasta;opts...)) for _ in 1:repeats]
    fwd = findrna(fasta;opts...,scoring=:forward)        # nhmmer-style rescoring
    ft = [@elapsed(findrna(fasta;opts...,scoring=:forward)) for _ in 1:repeats]
    nt = [@elapsed(external(nhcmd,joinpath(out,"nhmmer.txt"),joinpath(out,"nhmmer.stderr"))) for _ in 1:repeats]
    bt = [@elapsed(external(bacmd,gff,joinpath(out,"barrnap.stderr"))) for _ in 1:repeats]
    write_gff(joinpath(out,"findrna.gff"),hits)
    refs = parse_nhmmer(tbl); br = parse_gff(gff)
    write_differences(joinpath(out,"nhmmer-differences.tsv"),hits,refs)
    write_differences(joinpath(out,"barrnap-differences.tsv"),hits,br)
    report = Dict("julia"=>string(VERSION),"threads"=>threads,"repeats"=>repeats,
        "cpu"=>Sys.CPU_NAME,"machine"=>Sys.MACHINE,
        "engine_source_sha256"=>bytes2hex(sha256(join(read.(sort(readdir(dirname(pathof(ProfileHMMs));join=true)),String)))),
        "input_sha256"=>bytes2hex(sha256(read(fasta))),"model_sha256"=>bytes2hex(sha256(read(hmm))),
        "findrna_first_call_seconds"=>cold,"findrna_warm_seconds"=>jt,"findrna_forward_warm_seconds"=>ft,
        "nhmmer_process_seconds"=>nt,"barrnap_process_seconds"=>bt,
        "findrna_hits"=>length(hits),"nhmmer_hits"=>length(refs),"barrnap_hits"=>length(br))
    open(joinpath(out,"timings.toml"),"w") do io; TOML.print(io,report); end
    write_scores(joinpath(out,"nhmmer-scores.tsv"),fwd,refs)
    println("Median seconds: RibosomalRNAFinder warm=",median_time(jt),"; forward=",median_time(ft),"; nhmmer=",median_time(nt),"; barrnap=",median_time(bt))
    if get(ENV,"BARRNAP_MODERN","false") == "true"
        # Different models/algorithms: report all raw calls; no fabricated parity.
        @testset "HMMER reference" begin
            @test !isempty(refs)
            compare_hits(hits,refs)
            compare_hits(fwd,refs)
            compare_scores(fwd,refs)
        end
        if isempty(ARGS)
            @testset "Modern Barrnap synthetic fixture" begin
                @test !isempty(br)
                compare_hits(hits,br)
            end
        else
            println("Modern Barrnap real-genome calls are reported separately; CM/HMM parity is not asserted.")
        end
    else
        occursin("0.9",bv) || error("Expected barrnap 0.9; set BARRNAP_MODERN=true for CM comparison")
        @testset "Reference hit identity, strand and endpoints" begin
            @test !isempty(refs)
            @test !isempty(br)
            compare_hits(hits,refs)
            compare_hits(hits,br)
            compare_hits(fwd,refs)
            compare_scores(fwd,refs)
        end
    end
end
abspath(PROGRAM_FILE) == (@__FILE__) && benchmark_main()
