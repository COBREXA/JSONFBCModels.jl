
# JSONFBCModels.jl -- Import and export of JSON flux-balanced constrained models


| Build status |
|:---:|
| [![CI](https://github.com/COBREXA/JSONFBCModels.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/COBREXA/JSONFBCModels.jl/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/COBREXA/JSONFBCModels.jl/branch/master/graph/badge.svg?token=A2ui7exGIH)](https://codecov.io/gh/COBREXA/JSONFBCModels.jl) |

This package defines an instance of the `AbstractFBCModel` interface (from the
package [AbstractFBCModels.jl](https://github.com/COBREXA/AbstractFBCModels.jl))
for models saved as JSON files (typically produced by
[cobrapy](https://opencobra.github.io/cobrapy/)). This allows you to easily use
the JSON formatted models in constraint-based modeling packages, and convert
them to other constraint-based metabolic modeling data formats.

The primary purpose of this is to provide JSON loading functionality for
[COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl) and
[FBCModelTests.jl](https://github.com/LCSB-BioCore/FBCModelTests.jl), but is
otherwise completely generic and can be used independently of these packages.

You should be able to load JSON models (called `JSONFBCModel`) via the
AbstractFBCModels interface:

```julia
import AbstractFBCModels as A
import JSONFBCModels

model = A.load("my_model.json")
A.reactions(model) # accessors defined by the interface (NB: access through the overloads: A.xyz)
```
See the docstrings of all the accessors available by entering `A.accessors()` in
the REPL. For example:

#### Acknowledgements

`JSONFBCModels.jl` was developed at the Luxembourg Centre for Systems
Biomedicine of the University of Luxembourg
([uni.lu/lcsb](https://www.uni.lu/lcsb))
and at Institute for Quantitative and Theoretical Biology at Heinrich Heine
University Düsseldorf ([qtb.hhu.de](https://www.qtb.hhu.de/en/)).
The development was supported by European Union's Horizon 2020 Programme under
PerMedCoE project ([permedcoe.eu](https://www.permedcoe.eu/)),
agreement no. 951773.

<img src="docs/src/assets/unilu.svg" alt="Uni.lu logo" height="64px">   <img src="docs/src/assets/lcsb.svg" alt="LCSB logo" height="64px">   <img src="docs/src/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/permedcoe.svg" alt="PerMedCoE logo" height="64px">
