JSONFBCModel(json::Dict{String,Any}) = begin
    rkey = first(intersect(keys(json), constants.keynames.reactions))
    isnothing(rkey) && throw(DomainError(keys(json), "JSON model has no reaction keys"))
    rs = json[rkey]

    mkey = first(intersect(keys(json), constants.keynames.metabolites))
    ms = json[mkey]
    isnothing(mkey) &&
        throw(DomainError(keys(json), "JSON model has no metabolite keys"))

    gkey = first(intersect(keys(json), constants.keynames.genes))
    gs = isnothing(gkey) ? [] : json[gkey]

    JSONFBCModel(
        json,
        Dict(extract_json_reaction_id(r, i) => i for (i, r) in enumerate(rs)),
        rs,
        Dict(extract_json_metabolite_id(m, i) => i for (i, m) in enumerate(ms)),
        ms,
        Dict(extract_json_gene_id(g, i) => i for (i, g) in enumerate(gs)),
        gs,
    )
end

# model
A.reactions(model::JSONFBCModel) =
    String[extract_json_reaction_id(r, i) for (i, r) in enumerate(model.reactions)]

A.n_reactions(model::JSONFBCModel) = length(model.reactions)

A.metabolites(model::JSONFBCModel) =
    String[extract_json_metabolite_id(m, i) for (i, m) in enumerate(model.metabolites)]

A.n_metabolites(model::JSONFBCModel) = length(model.metabolites)

A.genes(model::JSONFBCModel) =
    String[extract_json_gene_id(g, i) for (i, g) in enumerate(model.genes)]

A.n_genes(model::JSONFBCModel) = length(model.genes)

function A.stoichiometry(model::JSONFBCModel)
    rxn_ids = A.reactions(model)
    met_ids = A.metabolites(model)

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
    Float64[
        get(rxn, "lower_bound", -constants.default_reaction_bound) for
        rxn in model.reactions
    ],
    Float64[
        get(rxn, "upper_bound", constants.default_reaction_bound) for rxn in model.reactions
    ],
)

A.balance(model::JSONFBCModel) = sparsevec(spzeros(A.n_metabolites(model)))

A.objective(model::JSONFBCModel) = sparsevec(
    Float64[float(get(rxn, "objective_coefficient", 0.0)) for rxn in model.reactions],
)

function A.reaction_gene_products_available(
    model::JSONFBCModel,
    rid::String,
    available::Function,
)
    x = get(model.reactions[model.reaction_index[rid]], "gene_reaction_rule", nothing)
    isnothing(x) && return nothing
    x = parse_gene_association(x)
    isnothing(x) && return nothing
    eval_gene_association(x, available)
end

function A.reaction_gene_association_dnf(model::JSONFBCModel, rid::String)
    x = get(model.reactions[model.reaction_index[rid]], "gene_reaction_rule", nothing)
    isnothing(x) && return nothing
    x = parse_gene_association(x)
    isnothing(x) && return nothing
    flatten_gene_association(x)
end

A.metabolite_formula(model::JSONFBCModel, mid::String) =
    parse_formula(get(model.metabolites[model.metabolite_index[mid]], "formula", nothing))

A.metabolite_charge(model::JSONFBCModel, mid::String) =
    parse_charge(get(model.metabolites[model.metabolite_index[mid]], "charge", nothing))

A.gene_annotations(model::JSONFBCModel, gid::String) = parse_annotations_or_notes(
    get(model.genes[model.gene_index[gid]], "annotation", A.Annotations()),
)

A.gene_notes(model::JSONFBCModel, gid::String) =
    parse_annotations_or_notes(get(model.genes[model.gene_index[gid]], "notes", A.Notes()))

A.reaction_annotations(model::JSONFBCModel, rid::String) = parse_annotations_or_notes(
    get(model.reactions[model.reaction_index[rid]], "annotation", A.Annotations()),
)

A.reaction_notes(model::JSONFBCModel, rid::String) = parse_annotations_or_notes(
    get(model.reactions[model.reaction_index[rid]], "notes", A.Notes()),
)

A.metabolite_annotations(model::JSONFBCModel, mid::String) = parse_annotations_or_notes(
    get(model.metabolites[model.metabolite_index[mid]], "annotation", A.Annotations()),
)

A.metabolite_notes(model::JSONFBCModel, mid::String) = parse_annotations_or_notes(
    get(model.metabolites[model.metabolite_index[mid]], "notes", A.Notes()),
)

A.reaction_name(model::JSONFBCModel, rid::String) =
    parse_name(get(model.reactions[model.reaction_index[rid]], "name", nothing))

A.metabolite_name(model::JSONFBCModel, mid::String) =
    parse_name(get(model.metabolites[model.metabolite_index[mid]], "name", nothing))

A.gene_name(model::JSONFBCModel, gid::String) =
    parse_name(get(model.genes[model.gene_index[gid]], "name", nothing))

A.reaction_stoichiometry(model::JSONFBCModel, rid::String)::Dict{String,Float64} =
    model.reactions[model.reaction_index[rid]]["metabolites"]

function A.metabolite_compartment(model::JSONFBCModel, mid::String)
    x = get(model.metabolites[model.metabolite_index[mid]], "compartment", nothing)
    return isa(x, String) ? x : nothing
end

function Base.convert(::Type{JSONFBCModel}, mm::A.AbstractFBCModel)
    if typeof(mm) == JSONFBCModel
        return mm
    end

    rxn_ids = A.reactions(mm)
    met_ids = A.metabolites(mm)
    gene_ids = A.genes(mm)
    S = A.stoichiometry(mm)
    lbs, ubs = A.bounds(mm)
    ocs = A.objective(mm)

    json = Dict{String,Any}()

    json[first(constants.keynames.genes)] = [
        Dict([
            "id" => gid,
            "name" => A.gene_name(mm, gid),
            "annotation" => A.gene_annotations(mm, gid),
            "notes" => A.gene_notes(mm, gid),
        ],) for gid in gene_ids
    ]

    json[first(constants.keynames.metabolites)] = [
        Dict([
            "id" => mid,
            "name" => A.metabolite_name(mm, mid),
            "formula" => unparse_formula(A.metabolite_formula(mm, mid)),
            "charge" => A.metabolite_charge(mm, mid),
            "compartment" => A.metabolite_compartment(mm, mid),
            "annotation" => A.metabolite_annotations(mm, mid),
            "notes" => A.metabolite_notes(mm, mid),
        ]) for mid in met_ids
    ]

    json[first(constants.keynames.reactions)] = [
        begin
            res = Dict{String,Any}()
            res["id"] = rid
            res["name"] = A.reaction_name(mm, rid)
            res["annotation"] = A.reaction_annotations(mm, rid)
            res["notes"] = A.reaction_notes(mm, rid)

            grr = A.reaction_gene_association_dnf(mm, rid)
            if !isnothing(grr)
                res["gene_reaction_rule"] = format_gene_association_dnf(grr)
            end

            res["lower_bound"] = lbs[ri]
            res["upper_bound"] = ubs[ri]
            res["objective_coefficient"] = ocs[ri]
            I, V = findnz(S[:, ri])
            res["metabolites"] =
                Dict{String,Float64}([met_ids[ii] => vv for (ii, vv) in zip(I, V)])
            identity(res)
        end for (ri, rid) in enumerate(rxn_ids)
    ]

    return JSONFBCModel(json)
end
