
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
