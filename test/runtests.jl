using Pkg

# Run each notebook in its own embedded Pluto environment, as Pluto does.
for name in (isempty(ARGS) ? ["amplitudes", "angles"] : ARGS)
    name in ("amplitudes", "angles") || error("Unknown notebook: $name")
    notebook = joinpath(@__DIR__, "..", "notebooks", "Jpsipipi_$(name).jl")
    source = read(notebook, String)
    mktempdir() do directory
        environment = realpath(directory)
        for (variable, file) in (
            ("PLUTO_PROJECT_TOML_CONTENTS", "Project.toml"),
            ("PLUTO_MANIFEST_TOML_CONTENTS", "Manifest.toml"),
        )
            contents = match(Regex(variable * " = \"\"\"(.*?)\"\"\"", "s"), source)
            isnothing(contents) && error("Missing $variable in $notebook")
            write(joinpath(environment, file), contents.captures[1])
        end
        Pkg.activate(environment)
        Pkg.instantiate()
        check = joinpath(@__DIR__, "check_notebook.jl")
        run(`$(Base.julia_cmd()) --project=$environment $check $name`)
    end
end
