@testset "IO" begin
    model = J.load(J.JSONFBCModel, iml1515_path)

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

    @test all(in.(J.filename_extensions(J.JSONFBCModel), Ref([
        "json"
        "JSON"
    ])))

    # TODO save
end
