@testset "IO" begin
    model = A.load(J.JSONFBCModel, iml1515_path)

    @test all(
        in.(
            keys(model.json),
            Ref([
                "metabolites"
                "id"
                "compartments"
                "reactions"
                "version"
                "genes"
            ]),
        ),
    )

    @test all(in.(A.filename_extensions(J.JSONFBCModel), Ref([
        "json"
        "JSON"
    ])))

    saved_path = joinpath(mktempdir(), "test-model.json")
    A.save(model, saved_path)
    test_iml1515_details(saved_path)   
end
