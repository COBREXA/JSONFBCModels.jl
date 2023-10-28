using Test
import AbstractFBCModels as A
import JSONFBCModels as J

include("test_utils.jl")

isdir("downloaded") || mkdir("downloaded")

iml1515_path = A.download_data_file(
    "http://bigg.ucsd.edu/static/models/iML1515.json",
    joinpath("downloaded", "iML1515.json"),
    "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
)

@testset "JSONFBCModels" begin
    @testset "Interface implemented" begin
        A.run_fbcmodel_type_tests(J.JSONFBCModel)
    end

    @testset "IO" begin
        include("io.jl")
        # TODO test convert

    end

    @testset "Accessors" begin
        A.run_fbcmodel_file_tests(J.JSONFBCModel, iml1515_path; name="iML1515", test_save=true)

        @testset "Specific details of iML1515" begin
            test_iml1515_details(iml1515_path)            
        end
    end
end
