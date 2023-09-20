using Test
import RequiredInterfaces as R
import AbstractFBCModels as A
import JSONFBCModels as J

isdir("downloaded") || mkdir("downloaded")

iml1515_path = A.download_data_file(
    "http://bigg.ucsd.edu/static/models/iML1515.json",
    joinpath("downloaded", "iML1515.json"),
    "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
)

@testset "JSONFBCModels" begin
    @testset "Interface implemented" R.check_implementations(A.AbstractFBCModel)
    # IO
    include("io.jl")
    # accessors
    include("accessors.jl")
end
