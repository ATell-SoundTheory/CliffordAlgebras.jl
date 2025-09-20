using Test

if VERSION >= v"1.8"
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
    @info "Skipping Aqua tests on Julia < 1.8"
end
