
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
