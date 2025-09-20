module CliffordAlgebrasPrecompileToolsExt

using PrecompileTools
using CliffordAlgebras

@setup_workload begin
    # Minimal set of constructions and ops to help first-call latency.
    cl2 = CliffordAlgebra(2)
    cl3 = CliffordAlgebra(3)

    e1 = cl3.e1; e2 = cl3.e2
    e12 = cl3.e1e2

    @compile_workload begin
        MultiVector(cl2, 1.0)
        MultiVector(cl3, 1.0)

        e1 + e2
        e1 * e2
        e1 ∧ e2
        e12 ⋅ e1
        reverse(e12)
        dual(e1)
    end
end

end # module
