
_json_rxn_name(r, i) = string(get(r, "id", "rxn$i"))

_json_met_name(m, i) = string(get(m, "id", "met$i"))

_json_gene_name(g, i) = string(get(g, "id", "gene$i"))

function parse_grr(str::Maybe{String})
    isnothing(str) && return nothing
    isempty("") && return nothing
    
    dnf = A.GeneAssociationDNF()
    for isozyme in string.(split(str, " or "))
        push!(
            dnf,
            string.(split(replace(isozyme, "(" => "", ")" => "", " and " => " "), " ")),
        )
    end
    return dnf
end

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

function parse_charge(x)
    if isa(x, Int)
        return x
    elseif isa(x, Float64)
        return Int(x)::Int
    elseif isa(x, String)
        return parse(Int, x)
    else
        nothing
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
    join([k * string(x[k]) for k in ks])
end

function unparse_grr(xs::Maybe{A.GeneAssociationDNF})
    isnothing(xs) && return nothing
    join([join(x, " and ") for x in xs], " or ")
end
