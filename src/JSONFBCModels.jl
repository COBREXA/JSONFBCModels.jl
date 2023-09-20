module JSONFBCModels

using DocStringExtensions
import AbstractFBCModels as A
using JSON, SparseArrays, Reexport

module Internal
using DocStringExtensions
import ..A
include("constants.jl")
include("utils.jl")
end

import .Internal

include("model.jl")
include("io.jl")

export JSONFBCModel

@reexport using AbstractFBCModels

end
