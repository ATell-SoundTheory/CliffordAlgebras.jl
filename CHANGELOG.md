# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog and this project adheres to Semantic Versioning.

## 0.1.4 - 2025-09-21
- Feature: PrettyTables integration via optional package extension; core now has a Unicode fallback renderer for `signaturetable` and `cayleytable`.
- Packaging: Move PrettyTables to `[weakdeps]` with `[extensions]` mapping; keep compat and add to `[extras]` to satisfy Pkg 1.8 validation.
- Docs: Mention PrettyTables extension usage and fallback in README and Algebra Catalog.
- CI: Fix Aqua 1.8 job by adjusting compat/extras; ensure Runtests, Documenter, Aqua, BenchSmoke all pass.

## 0.1.3 - 2025-09-20
- Performance: cache multiplication tables per algebra type (reduces repeated-product overhead)
- API: add `hash` and `isequal` for `MultiVector` (dict/set-friendly, stable hashing)
- Quality: integrate Aqua.jl; fix unbound type parameter warnings in constructors
- Precompile: add optional PrecompileTools extension (keeps Julia 1.6 CI green)
- Compat: add `[compat] PrecompileTools = "1"` (weakdep); relax PrettyTables = 2, StaticArrays = 1
- Docs: enable doctests; improve tutorial/index; add Developer Guide
- CI: modernize actions (checkout v4); test on 1.6, 1.10, latest 1.x; Documenter build fixes

## 0.1.2
- Current published version
