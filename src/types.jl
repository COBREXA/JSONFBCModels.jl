
"""
$(TYPEDEF)

A representation of the contents of a JSON model, i.e. a metabolic model read
from a file ending with `.json`. The data is stored in the Julia version of
JSON as given by `JSON.jl`, basically as vectors and dictionaries.

Because the JSON format is very open, some of the data in the JSON models are
not represented by the AbstractFBCModel interface, and thus gets lost upon
conversion. Always check the conversion results to avoid losing information.

Since direct work with the JSON structure is not very efficient, the model
structure caches some of the internal structure in the extra fields.  The
single-parameter [`JSONFBCModel`](@ref) constructor creates these caches
automatically from the input JSON. The model structure is designed as
read-only, and any changes in the JSON fields usually break the model.

# Fields
$(TYPEDFIELDS)
"""
struct JSONFBCModel <: A.AbstractFBCModel
    json::JSON.Object{String,Any}
    reaction_index::Dict{String,Int}
    reactions::Vector{Any}
    metabolite_index::Dict{String,Int}
    metabolites::Vector{Any}
    gene_index::Dict{String,Int}
    genes::Vector{Any}
end

const Maybe = A.Maybe
