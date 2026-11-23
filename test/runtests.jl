using Test, RibosomalRNAFinder, ProfileHMMs, BioSequences, Tables
include("fixtures.jl")

@testset "RibosomalRNAFinder" begin
    ps = parse_hmm_file(default_profiles())
    @test getfield.(ps, :name) == ["16S_rRNA", "23S_rRNA", "5S_rRNA"]
    @test getfield.(ps, :length) == [1533, 2906, 119]
    mktempdir() do dir
        path = joinpath(dir, "fixture.fna")
        _, seqs = make_fixture(path)
        r = findrna(path; progress=false)
        @test r isa RNAResult
        @test length(r) == 9
        @test all(h -> h.seqid ∉ ("negative", "ambiguous"), r)
        @test !any(h -> h.partial, r)
        for p in ps, (suffix, strand) in (("forward", '+'), ("reverse", '-'))
            hits = filter(h -> h.seqid == p.name * "_" * suffix && h.model == p.name, r)
            @test length(hits) == 1
            if length(hits) == 1
                @test hits[1].strand == strand
                @test abs(hits[1].start - 81) <= 5
                @test abs(hits[1].stop - (80 + p.length)) <= 5
            end
        end
        @test length(filter(h -> h.seqid == "two_16S", r)) == 2
        real5s = filter(h -> h.seqid == "real_5S", r)
        @test [(h.model, h.start, h.stop, h.strand) for h in real5s] == [("5S_rRNA", 2401, 2510, '+')]

        # Tables.jl, indexing, extraction and identifiers
        @test Tables.istable(typeof(r))
        cols = Tables.columntable(r)
        @test cols.seqid == [h.seqid for h in r] && eltype(cols.partial) == Bool
        @test r[1:2] isa RNAResult && length(r[1:2]) == 2
        @test filter(h -> h.model == "5S_rRNA", r) isa RNAResult
        extracted = extract_sequences(r)
        @test length(extracted) == length(r)
        @test all(length(s) == length(h) for (s, h) in zip(extracted, r))
        rev = only(filter(h -> h.seqid == "16S_rRNA_reverse", r))
        fwd = only(filter(h -> h.seqid == "16S_rRNA_forward", r))
        @test extract_sequences(r[[findfirst(==(rev), r)]])[1] ==
              extract_sequences(r[[findfirst(==(fwd), r)]])[1]
        ids = hit_ids(r)
        @test length(unique(ids)) == length(r)
        @test ids[1] == "$(r[1].seqid):$(r[1].start)-$(r[1].stop)($(r[1].strand))"

        # Display and file output
        shown = sprint(show, MIME"text/plain"(), r)
        @test !occursin('\e', shown) && occursin("RNAResult", shown) && occursin("16S_rRNA × 4", shown)
        gff = read(write_gff(joinpath(dir, "hits.gff"), r), String)
        @test startswith(gff, "##gff-version 3\n")
        @test occursin("Name=16S_rRNA;product=16S ribosomal RNA", gff)
        @test count(==('\n'), gff) == length(r) + 1
        fa = readlines(write_fasta(joinpath(dir, "hits.fa"), r))
        @test count(startswith(">"), fa) == length(r)
        @test fa[1] == ">" * ids[1] * " " * r[1].model

        # Options
        only16 = findrna(path; models=["16S_rRNA"], progress=false)
        @test all(h -> h.model == "16S_rRNA", only16) && length(only16) == 4
        @test_throws ArgumentError findrna(path; models=["no_such_model"], progress=false)
        @test_throws ArgumentError findrna(path; models="16S_rRNA", progress=false)
        @test_throws ArgumentError findrna(path; reject=2, progress=false)
        @test_throws ArgumentError findrna(joinpath(dir, "absent.fna"); progress=false)
        empty = findrna(path; models=String[], progress=false)
        @test isempty(empty) && isempty(extract_sequences(empty))
        strict = findrna(path; evalue=1e-40, progress=false)   # synthetic 5S: E ≈ 1e-34
        @test !isempty(strict) && all(h -> h.model != "5S_rRNA", strict)
        @test RibosomalRNAFinder._gff_seqid("contig 1;x") == "contig%201%3Bx"
        @test RibosomalRNAFinder._gff_attr("a b;c=d,é") == "a b%3Bc%3Dd%2Cé"
    end
end

@testset "Partial and rejected hits (barrnap rules)" begin
    p16 = parse_hmm_file(default_profiles())[1]
    cons = LongDNA{4}(join("ACGT"[argmin(p16.match_emissions[:, k])] for k in 1:p16.length))
    mktempdir() do dir
        path = joinpath(dir, "fragments.fna")
        open(path, "w") do io
            println(io, ">half\n", "N"^50, cons[1:760], "N"^50)      # 50% of the model: partial
            println(io, ">tenth\n", "N"^50, cons[1:150], "N"^50)     # 10%: rejected
        end
        r = findrna(path; models=["16S_rRNA"], progress=false)
        @test [h.seqid for h in r] == ["half"]
        @test only(r).partial
        @test occursin("(partial)", sprint(write_gff, r))
        kept = findrna(path; models=["16S_rRNA"], reject=0.0, progress=false)
        @test sort([h.seqid for h in kept]) == ["half", "tenth"]
        @test !findrna(path; models=["16S_rRNA"], lencutoff=0.4, progress=false)[1].partial
    end
end

@testset "nhmmer-style scoring (scoring=:forward)" begin
    mktempdir() do dir
        path = joinpath(dir, "fixture.fna")
        make_fixture(path)
        vit = findrna(path; progress=false)
        fwd = findrna(path; progress=false, scoring=:forward)
        key(h) = (h.seqid, h.model, h.start, h.stop, h.strand)
        @test key.(fwd) == key.(vit)
        @test all(h -> h.bias == 0, vit) && all(h -> h.bias >= 0, fwd)
        # nhmmer 3.4 scores this 5S (CP148752.1:116377-116486) at 51.6 bits.
        real5s = only(filter(h -> h.seqid == "real_5S", fwd))
        @test abs(real5s.bits - 51.6) <= 2
        @test_throws ArgumentError findrna(path; progress=false, scoring=:fast)
    end
end
