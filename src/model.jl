"""
$(TYPEDEF)

A struct used to store the contents of a JSON model, i.e. a model read from a
file ending with `.json`. These model files typically store all the model data
in arrays of JSON objects (represented in Julia as vectors of dictionaries).

Usually, not all of the fields of the input JSON can be easily represented when
converting to other models, care should be taken to avoid losing information.

Direct work with the `json` structure is not very efficient; the model
structure therefore caches some of the internal structure in the extra fields.
The single-parameter [`JSONFBCModel`](@ref) constructor creates these caches
correctly from the `json`. The model structure is designed as read-only, and
changes in `json` invalidate the cache.

# Example
````
model = load("some_model.json")
model.json # see the actual underlying JSON
````

# Fields
$(TYPEDFIELDS)
"""
struct JSONFBCModel <: A.AbstractFBCModel
    json::Dict{String,Any}
    reaction_index::Dict{String,Int}
    reactions::Vector{Any}
    metabolite_index::Dict{String,Int}
    metabolites::Vector{Any}
    gene_index::Dict{String,Int}
    genes::Vector{Any}
end

JSONFBCModel(json::Dict{String,Any}) = begin
    rkey = A.guesskey(keys(json), Internal.constants.keynames.reactions)
    isnothing(rkey) && throw(DomainError(keys(json), "JSON model has no reaction keys"))
    rs = json[rkey]

    mkey = A.guesskey(keys(json), Internal.constants.keynames.metabolites)
    ms = json[mkey]
    isnothing(mkey) && throw(DomainError(keys(json), "JSON model has no metabolite keys"))

    gkey = A.guesskey(keys(json), Internal.constants.keynames.genes)
    gs = isnothing(gkey) ? [] : json[gkey]

    JSONFBCModel(
        json,
        Dict(Internal._json_rxn_name(r, i) => i for (i, r) in enumerate(rs)),
        rs,
        Dict(Internal._json_met_name(m, i) => i for (i, m) in enumerate(ms)),
        ms,
        Dict(Internal._json_gene_name(g, i) => i for (i, g) in enumerate(gs)),
        gs,
    )
end

A.reactions(model::JSONFBCModel) =
    [Internal._json_rxn_name(r, i) for (i, r) in enumerate(model.reactions)]

A.n_reactions(model::JSONFBCModel) = length(model.reactions)

A.metabolites(model::JSONFBCModel) =
    [Internal._json_met_name(m, i) for (i, m) in enumerate(model.metabolites)]

A.n_metabolites(model::JSONFBCModel) = length(model.metabolites)

A.genes(model::JSONFBCModel) =
    [Internal._json_gene_name(g, i) for (i, g) in enumerate(model.genes)]

A.n_genes(model::JSONFBCModel) = length(model.genes)

function A.stoichiometry(model::JSONFBCModel)
    rxn_ids = reactions(model)
    met_ids = metabolites(model)

    n_entries = 0
    for r in model.reactions
        for _ in r["metabolites"]
            n_entries += 1
        end
    end

    MI = Vector{Int}()
    RI = Vector{Int}()
    SV = Vector{Float64}()
    sizehint!(MI, n_entries)
    sizehint!(RI, n_entries)
    sizehint!(SV, n_entries)

    for (i, rid) in enumerate(rxn_ids)
        r = model.reactions[model.reaction_index[rid]]
        for (mid, coeff) in r["metabolites"]
            haskey(model.metabolite_index, mid) || throw(
                DomainError(
                    met_id,
                    "Unknown metabolite found in stoichiometry of $(rxn_ids[i])",
                ),
            )

            push!(MI, model.metabolite_index[mid])
            push!(RI, i)
            push!(SV, coeff)
        end
    end
    return SparseArrays.sparse(MI, RI, SV, length(met_ids), length(rxn_ids))
end

A.bounds(model::JSONFBCModel) = (
    [
        get(rxn, "lower_bound", -Internal.constants.default_reaction_bound) for
        rxn in model.reactions
    ],
    [
        get(rxn, "upper_bound", Internal.constants.default_reaction_bound) for
        rxn in model.reactions
    ],
)

A.balance(model::JSONFBCModel) = spzeros(n_metabolites(model))

A.objective(model::JSONFBCModel) =
    sparse([float(get(rxn, "objective_coefficient", 0.0)) for rxn in model.reactions])

A.reaction_gene_associations(model::JSONFBCModel, rid::String) = A.maybemap(
    A.parse_grr,
    get(model.reactions[model.reaction_index[rid]], "gene_reaction_rule", nothing),
)

A.metabolite_formula(model::JSONFBCModel, mid::String) = A.maybemap(
    A.parse_formula,
    get(model.metabolites[model.metabolite_index[mid]], "formula", nothing),
)

A.metabolite_charge(model::JSONFBCModel, mid::String) =
    get(model.metabolites[model.metabolite_index[mid]], "charge", 0)

A.metabolite_compartment(model::JSONFBCModel, mid::String) =
    get(model.metabolites[model.metabolite_index[mid]], "compartment", nothing)

A.gene_annotations(model::JSONFBCModel, gid::String) = A.maybemap(
    Internal._parse_annotations,
    get(model.genes[model.gene_index[gid]], "annotation", nothing),
)

A.gene_notes(model::JSONFBCModel, gid::String) = A.maybemap(
    Internal._parse_notes,
    get(model.genes[model.gene_index[gid]], "notes", nothing),
)

A.reaction_annotations(model::JSONFBCModel, rid::String) = A.maybemap(
    Internal._parse_annotations,
    get(model.reactions[model.reaction_index[rid]], "annotation", nothing),
)

A.reaction_notes(model::JSONFBCModel, rid::String) = A.maybemap(
    Internal._parse_notes,
    get(model.reactions[model.reaction_index[rid]], "notes", nothing),
)

A.metabolite_annotations(model::JSONFBCModel, mid::String) = A.maybemap(
    Internal._parse_annotations,
    get(model.metabolites[model.metabolite_index[mid]], "annotation", nothing),
)

A.metabolite_notes(model::JSONFBCModel, mid::String) = A.maybemap(
    Internal._parse_notes,
    get(model.metabolites[model.metabolite_index[mid]], "notes", nothing),
)

A.reaction_stoichiometry(model::JSONFBCModel, rid::String) =
    model.reactions[model.reaction_index[rid]]["metabolites"]

A.reaction_name(model::JSONFBCModel, rid::String) =
    get(model.reactions[model.reaction_index[rid]], "name", nothing)

A.metabolite_name(model::JSONFBCModel, mid::String) =
    get(model.metabolites[model.metabolite_index[mid]], "name", nothing)

A.gene_name(model::JSONFBCModel, gid::String) =
    get(model.genes[model.gene_index[gid]], "name", nothing)

function Base.convert(::Type{JSONFBCModel}, mm::AbstractFBCModel)
    if typeof(mm) == JSONFBCModel
        return mm
    end

    rxn_ids = variable_ids(mm)
    met_ids = metabolite_ids(mm)
    gene_ids = genes(mm)
    S = stoichiometry(mm)
    lbs, ubs = variable_bounds(mm)
    ocs = objective(mm)

    json = Dict{String,Any}()

    json["annotation"] = model_annotations(mm)
    json["notes"] = model_notes(mm)

    json[first(constants.keynames.genes)] = [
        Dict([
            "id" => gid,
            "name" => gene_name(mm, gid),
            "annotation" => gene_annotations(mm, gid),
            "notes" => gene_notes(mm, gid),
        ],) for gid in gene_ids
    ]

    json[first(constants.keynames.mets)] = [
        Dict([
            "id" => mid,
            "name" => metabolite_name(mm, mid),
            "formula" => maybemap(unparse_formula, metabolite_formula(mm, mid)),
            "charge" => metabolite_charge(mm, mid),
            "compartment" => metabolite_compartment(mm, mid),
            "annotation" => metabolite_annotations(mm, mid),
            "notes" => metabolite_notes(mm, mid),
        ]) for mid in met_ids
    ]

    json[first(constants.keynames.rxns)] = [
        begin
            res = Dict{String,Any}()
            res["id"] = rid
            res["name"] = reaction_name(mm, rid)
            res["subsystem"] = reaction_subsystem(mm, rid)
            res["annotation"] = reaction_annotations(mm, rid)
            res["notes"] = reaction_notes(mm, rid)

            grr = reaction_gene_associations(mm, rid)
            if !isnothing(grr)
                res["gene_reaction_rule"] = unparse_grr(String, grr)
            end

            res["lower_bound"] = lbs[ri]
            res["upper_bound"] = ubs[ri]
            res["objective_coefficient"] = ocs[ri]
            I, V = findnz(S[:, ri])
            res["metabolites"] =
                Dict{String,Float64}([met_ids[ii] => vv for (ii, vv) in zip(I, V)])
            res
        end for (ri, rid) in enumerate(rxn_ids)
    ]

    return JSONFBCModel(json)
end
