using Documenter
using JuliQAOA

using Graphs

makedocs(
    modules = [JuliQAOA],
    sitename = "JuliQAOA",
    authors = "John Golden",
    pages = [
        "Home" => "index.md",
        "Basic Use" => "basics.md",
        "Manual" => ["mixers.md",
                     "cost_funcs.md",
                     "eval.md",
                     "angle_finding.md",
                     "utils.md"],
        "Examples" => "examples.md"
    ],
    strict = false
)

deploydocs(
    repo = "github.com/lanl/JuliQAOA.jl.git",
)