
_json_rxn_name(r, i) = string(get(r, "id", "rxn$i"))

_json_met_name(m, i) = string(get(m, "id", "met$i"))

_json_gene_name(g, i) = string(get(g, "id", "gene$i"))

function _parse_annotations(x)
    A.Annotations([k => if typeof(v) == String
        [v]
    else
        convert(Vector{String}, v)
    end for (k, v) in x])
end

_parse_notes(x) = _parse_annotations(x)

"""
$(TYPEDSIGNATURES)

Unfortunately, many model types that contain dictionares do not have
standardized field names, so we need to try a few possibilities and guess the
best one. The keys used to look for valid field names should be ideally
specified as constants in `src/base/constants.jl`.
"""
function guesskey(avail, possibilities)
    x = intersect(possibilities, avail)

    if isempty(x)
        @debug "could not find any of keys: $possibilities"
        return nothing
    end

    if length(x) > 1
        @debug "Possible ambiguity between keys: $x"
    end
    return x[1]
end
