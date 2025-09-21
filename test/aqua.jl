using Test

# Run Aqua locally on Julia >= 1.8. In CI, only run Aqua if AQUA_ONLY=true to avoid
# duplicating checks in the main test matrix. This keeps 1.6 CI clean and Aqua isolated.
if VERSION >= v"1.8" && (get(ENV, "CI", "") != "true" || get(ENV, "AQUA_ONLY", "false") == "true")
    using Aqua
    using CliffordAlgebras
    @testset "Aqua.jl" begin
        # Configure checks:
        # - unbound_args: outer constructors for parametric types can trigger false positives; mark as broken for now.
        # - deps_compat: skip extras check to avoid requiring compat for stdlib Test.
        Aqua.test_all(
            CliffordAlgebras;
            piracies=true,
            unbound_args=(broken=false,),
            deps_compat=(check_extras=false,),
        )
    end
else
    @info "Skipping Aqua tests (either Julia < 1.8 or CI main job without AQUA_ONLY)"
end
