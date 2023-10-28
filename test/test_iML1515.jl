
@testset "Test the expected contents of iML1515.json" begin

    # the path is inherited from the main testset
    model = A.load(JSONFBCModel, joinpath(@__DIR__, "test-models", "iML1515.json"))

    @test "SHK3Dr" in A.reactions(model)
    @test A.n_reactions(model) == 2712
    @test "co2_e" in A.metabolites(model)
    @test A.n_metabolites(model) == 1877
    @test "b0008" in A.genes(model)
    @test A.n_genes(model) == 1516
    @test all(size(A.stoichiometry(model)) .== (1877, 2712))
    @test all(length.(A.bounds(model)) .== 2712)
    @test all(A.balance(model) .== 0)
    @test A.objective(model)[2669] == 1
    @test all(
        in.(A.reaction_gene_association_dnf(model, "FBA"), Ref([["b2925"], ["b2097"]])),
    )
    @test A.metabolite_formula(model, "atp_c")["C"] == 10
    @test A.metabolite_charge(model, "atp_c") == -4
    @test A.metabolite_compartment(model, "atp_c") == "c"
    @test haskey(A.gene_annotations(model, "b0008"), "uniprot")
    @test haskey(A.gene_notes(model, "b0008"), "original_bigg_ids")
    @test haskey(A.reaction_annotations(model, "FBA"), "rhea")
    @test haskey(A.reaction_notes(model, "FBA"), "original_bigg_ids")
    @test haskey(A.metabolite_annotations(model, "atp_c"), "biocyc")
    @test haskey(A.metabolite_notes(model, "atp_c"), "original_bigg_ids")
    @test A.reaction_stoichiometry(model, "FBA")["dhap_c"] == 1.0
    @test A.reaction_name(model, "FBA") == "Fructose-bisphosphate aldolase"
    @test A.metabolite_name(model, "atp_c") == "ATP C10H12N5O13P3"
    @test A.gene_name(model, "b0008") == "talB"
end
