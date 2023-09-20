@testset "Accessors" begin
    model = J.load(J.JSONFBCModel, iml1515_path)

    @test "SHK3Dr" in J.reactions(model)
    @test J.n_reactions(model) == 2712
    @test "co2_e" in J.metabolites(model)
    @test J.n_metabolites(model) == 1877
    @test "b0008" in J.genes(model)
    @test J.n_genes(model) == 1516
    @test all(size(J.stoichiometry(model)) .== (1877, 2712))
    @test all(length.(J.bounds(model)) .== 2712)
    @test all(J.balance(model) .== 0)
    @test J.objective(model)[2669] == 1
    @test all(in.(J.reaction_gene_associations(model, "FBA"), Ref([["b2925"], ["b2097"]])))
    @test J.metabolite_formula(model, "atp_c")["C"] == 10
    @test J.metabolite_charge(model, "atp_c") == -4
    @test J.metabolite_compartment(model, "atp_c") == "c"
    @test haskey(J.gene_annotations(model, "b0008"), "uniprot")
    @test haskey(J.gene_notes(model, "b0008"), "original_bigg_ids")
    @test haskey(J.reaction_annotations(model, "FBA"), "rhea")
    @test haskey(J.reaction_notes(model, "FBA"), "original_bigg_ids")
    @test haskey(J.metabolite_annotations(model, "atp_c"), "biocyc")
    @test haskey(J.metabolite_notes(model, "atp_c"), "original_bigg_ids")
    @test J.reaction_stoichiometry(model, "FBA")["dhap_c"] == 1.0
    @test J.reaction_name(model, "FBA") == "Fructose-bisphosphate aldolase"
    @test J.metabolite_name(model, "atp_c") == "ATP C10H12N5O13P3"
    @test J.gene_name(model, "b0008") == "talB"
end
