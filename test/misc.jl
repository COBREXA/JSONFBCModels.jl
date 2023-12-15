
@testset "Miscellaneous tests" begin
    # the path is inherited from the main testset
    model = A.load(JSONFBCModel, joinpath(@__DIR__, "test-models", "iML1515.json"))
    # test that no superfluous keys are generated
    @test all(
        key in ["metabolites", "id", "compartments", "reactions", "version", "genes"] for
        key in keys(model.json)
    )

    # caps suffix is quite common
    @test all(suffix in A.filename_extensions(JSONFBCModel) for suffix in ["json", "JSON"])
end

@testset "Corner cases" begin
    import JSONFBCModels:
        eval_gene_association, flatten_gene_association, parse_charge, sortunique

    @test parse_charge(1) == 1
    @test parse_charge(2.0) == 2
    @test parse_charge("3") == 3
    @test parse_charge("4.0") == 4
    @test parse_charge(nothing) == nothing
    @test_throws ArgumentError parse_charge("totally positive charge")
    @test_throws DomainError parse_charge(["very charged"])
    @test_throws DomainError eval_gene_association(:(xor(gene("a"), gene("b"))), _ -> false)
    @test_throws DomainError flatten_gene_association(:(xor(gene("a"), gene("b"))))
    @test sortunique([3, 2, 2, 1]) == [1, 2, 3]
end
