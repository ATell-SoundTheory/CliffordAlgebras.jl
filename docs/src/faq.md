# FAQ & Tips

## How do I install the package?

This package is currently hosted on GitHub. Install with:

```julia
import Pkg
Pkg.add(url="https://github.com/ATell-SoundTheory/CliffordAlgebras.jl")
```

For local development:

```julia
import Pkg
Pkg.develop(path="/absolute/path/to/CliffordAlgebras.jl")
```

Minimum Julia version: 1.6.

## How do I type Unicode operators like ∧ ⋅ ⨼ ⨽ ×₋ ≀ ?

Use Julia's LaTeX-style tab completions in the REPL and many editors:

- `\wedge<Tab>` → `∧` (exterior)
- `\vee<Tab>` → `∨` (regressive)
- `\cdot<Tab>` → `⋅` (fat dot)
- `\lrcorner<Tab>` → `⨼` (left contraction)
- `\llcorner<Tab>` → `⨽` (right contraction)
- `\star<Tab>` → `⋆` (scalar)
- `\times<Tab>\_ -<Tab>` → `×₋` (commutator)
- `\times<Tab>\_ +<Tab>` → `×₊` (anti-commutator)
- `\Vert<Tab>` → `≀` (sandwich)

ASCII alternatives (function calls):

- `exteriorprod(a,b)` for `a ∧ b`
- `fatdotprod(a,b)` for `a ⋅ b`
- `leftcontractionprod(a,b)` for `a ⨼ b`
- `rightcontractionprod(a,b)` for `a ⨽ b`
- `scalarprod(a,b)` for `a ⋆ b`
- `commutatorprod(a,b)` for `a ×₋ b`
- `anticommutatorprod(a,b)` for `a ×₊ b`
- `sandwichproduct(a,b)` for `a ≀ b`

## What predefined algebras are available?

Common aliases include:

- `:Complex` (ℂ), `:Quaternions` (ℍ)
- `:Spacetime` (STA)
- `:PGA2D`, `:PGA3D`
- `:CGA2D`, `:CGA3D`

See the examples on the Home page for usage.

## How do I benchmark operations?

This repository includes a small benchmark suite under `benchmark/`.

- Quick run:

```julia
julia --project=benchmark benchmark/run.jl
```

- Full PkgBenchmark run:

```julia
julia --project=benchmark -e 'using PkgBenchmark, CliffordAlgebras; benchmarkpkg(CliffordAlgebras)'
```

## Can I use ASCII-only code?

Yes. All Unicode operators have function equivalents (see above). You can also import and call `CliffordAlgebras.<name>` functions directly.

## I get a MethodError mixing different algebras

Operations require both operands to belong to the same algebra. Construct both operands from the same `CliffordAlgebra(...)` instance.

## `prune` changes the type — is that intended?

Yes. `prune` removes zero coefficients and may return a different sparse type. Use `extend(mv)` if you need a full, fixed representation.
