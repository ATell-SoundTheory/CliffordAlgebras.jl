# Contributing to CliffordAlgebras.jl

Thanks for your interest in improving CliffordAlgebras.jl!

## Getting started

- Julia 1.6+ (1.10 and 1.latest are CI-tested)
- Clone the repo, then run tests:

```
julia --project -e 'using Pkg; Pkg.test()'
```

If you’re developing changes locally, it can be convenient to use a temp env and `Pkg.develop(path=...)` to isolate dependencies.

## Code style & performance

- Prefer type-stable APIs. Use `@inferred` in tests for hot paths.
- Leverage sparse multivector structure; avoid `extend` unless it’s necessary for stability.
- Be careful adding `@generated` functions—ensure they depend only on types, not runtime values.
- For new public behavior, add tests. For new algorithms, add a short note in the docs.

## Docs & doctests

- We enable doctests in Documenter. Prefer `jldoctest` blocks and smoke tests with stable outputs.
- Run docs locally:

```
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'
```

## CI

- The repo uses GitHub Actions for tests across Julia versions and OSes, and Documenter for docs.
- Please ensure CI is green before requesting review.

## Reporting issues

- Include minimal reproductions and the algebra signature in question (e.g., `CliffordAlgebra(3,0,0)`).
- For performance issues, include `@time`/`@btime` output and the sizes/grades involved.

## License

By contributing, you agree your contributions will be licensed under the project’s MIT license.
