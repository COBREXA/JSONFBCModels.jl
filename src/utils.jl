
extract_json_reaction_id(r, i) = string(get(r, "id", "rxn$i"))

extract_json_metabolite_id(m, i) = string(get(m, "id", "met$i"))

extract_json_gene_id(g, i) = string(get(g, "id", "gene$i"))

function parse_formula(x::Maybe{String})
    isnothing(x) && return nothing
    x == "" && return nothing

    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"
    for m in eachmatch(pattern, x)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end
    return res
end

function parse_charge(x)::Maybe{Int}
    if isa(x, Int)
        x
    elseif isa(x, Float64)
        Int(x)
    elseif isa(x, String)
        Int(parse(Float64, x))
    elseif isnothing(x)
        nothing
    else
        throw(DomainError(x, "cannot parse charge"))
    end
end

parse_name(x) = isa(x, String) ? x::String : nothing

function parse_annotations_or_notes(x)
    isa(x, A.Annotations) && return x

    a_or_n = A.Annotations()
    for (k, vs) in x
        if isa(vs, String)
            a_or_n[k] = String[vs]
        else
            a_or_n[k] = String[v for v in vs]
        end
    end
    return a_or_n
end

function unparse_formula(x::Maybe{A.MetaboliteFormula})
    isnothing(x) && return nothing
    ks = sort(collect(keys(x)))
    join(k * string(x[k]) for k in ks)
end
