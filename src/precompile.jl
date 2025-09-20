module __PrecompileWorkload

using PrecompileTools
using ..CliffordAlgebras

@setup_workload begin
    # Minimal set of constructions and ops to help first-call latency.
    # Keep this lightweight to avoid slowing down precompilation.
    cl2 = CliffordAlgebra(2)
    cl3 = CliffordAlgebra(3)
    pga = CliffordAlgebra(:PGA3D)

    e1 = cl3.e1; e2 = cl3.e2; e3 = cl3.e3
    e12 = cl3.e1e2; e23 = cl3.e2e3; e31 = cl3.e3e1

    @compile_workload begin
        # Constructors and basic ops
        MultiVector(cl2, 1.0)
        MultiVector(cl3, 1.0)
        MultiVector(pga, 1.0)

        e1 + e2
        e1 * e2
        e1 ∧ e2
    e12 ⋅ e3
    reverse(e12)
    dual(e1)
    end
end

end # module
