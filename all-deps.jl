using Pkg;
# ODE solver
Pkg.add("ODE")
# Python-based plotting, for more interactive and high-quality plots
Pkg.add("PyPlot")
# Adding Atom and VSCode deps
Pkg.develop("Revise")
Pkg.develop("Atom")
Pkg.develop("Juno")
Pkg.develop("DocumentFormat")
Pkg.develop("CSTParser")
Pkg.develop("StaticLint")
Pkg.develop("SymbolServer")
Pkg.develop("LanguageServer")
Pkg.develop("Linter")
