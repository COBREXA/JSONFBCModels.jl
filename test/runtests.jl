import AbstractFBCModels as A
import JSONFBCModels: JSONFBCModel

using Test

@testset "JSONFBCModels tests" begin
    A.run_fbcmodel_type_tests(JSONFBCModel)

    modeldir = joinpath(@__DIR__, "test-models")
    mkpath(modeldir)

    for (name, url, hash, ts) in [
        (
            "e_coli_core",
            "http://bigg.ucsd.edu/static/models/e_coli_core.json",
            "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
            true,
        ),
        (
            "iJO1366",
            "http://bigg.ucsd.edu/static/models/iJO1366.json",
            "9376a93f62ad430719f23e612154dd94c67e0d7c9545ed9d17a4d0c347672313",
            true,
        ),
        (
            "iML1515",
            "http://bigg.ucsd.edu/static/models/iML1515.json",
            "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
            true,
        ),
    ]
        path = joinpath(modeldir, "$name.json")
        A.download_data_file(url, path, hash)
        A.run_fbcmodel_file_tests(JSONFBCModel, path; name, test_save = ts)
    end

    include("test_iML1515.jl")
    include("misc.jl")
end
