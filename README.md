# Julia scripts

My Julia scripts used for numerical computation. I like [chaotic systems][1] and [Sturm-Liouville problems][2] the most, so they are the main ones I solve in it.

To (hopefully) install all dependencies required for the scripts in this repository, and for turning Atom and VSCode into suitable IDEs run [all-deps.jl](all-deps.jl) in Julia.

## [Airy Sturm-Liouville problem solved using rational transformation of Chebyshev extrema grid to [0,∞]][3]

View [airy-sle-problem-explanation.pdf][4] for a technical explanation of the problem and the method of numerically solving it.

[1]: https://en.wikipedia.org/wiki/Chaos_theory
[2]: https://en.wikipedia.org/wiki/Sturm%E2%80%93Liouville_theory
[3]: SLE/airy-rat.jl
[4]: SLE/airy-sle-problem-explanation.pdf
