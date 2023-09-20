
A.load(type::Type{JSONFBCModel}, filename::String) = JSONFBCModel(JSON.parsefile(filename))

A.save(model::JSONFBCModel, filename::String) =
    open(f -> JSON.print(f, model.json), filename, "w")

A.filename_extensions(type::Type{JSONFBCModel}) = ["json", "JSON"]
