using CliffordAlgebras
using Test
import LinearAlgebra.SingularException

@testset "CliffordAlgebras.jl" begin
    @testset "inference pga3d" begin
        pga = CliffordAlgebra(3,0,1)
        @inferred MultiVector(pga, 1.0)
        pt  = @inferred dual(pga.e4 + 3.2pga.e1 + 1.3pga.e2-4.3pga.e3)
        pt1 = @inferred dual(pga.e4 + pga.e1)
        pt2 = @inferred dual(pga.e4 + pga.e2 - pga.e3)
        l = @inferred ∨(pt1, pt2)
        @inferred norm(l)
        l = @inferred l / norm(l)
        @inferred CliffordAlgebras.exp_trig(-pi/6*l)
        motor1 = @inferred exp(-pi/6*l)
        motor2 = @inferred exp(-1/2*pga.e3e4)
        motor = @inferred motor1 * motor2

        @inferred Λᵏ(motor, Val(0))
        @inferred Λᵏ(motor, Val(1))
        @inferred Λᵏ(motor, Val(2))
        @inferred Λᵏ(motor, Val(3))
        @inferred Λᵏ(motor, Val(4))
        @inferred extend(motor)

        @inferred ≀(motor, pt)
        @inferred ≀(motor1, pt)
        @inferred ≀(motor2, pt)
        @inferred ≀(motor, l)
        @inferred ≀(motor, motor1)
    end
    @testset "isapprox" begin
        Cl = CliffordAlgebra
        @test Cl(1,1,1).e1 ≈ Cl(1,1,1).e1
        @test !(Cl(2,0,1).e1 ≈ Cl(1,1,1).e1)
        @test !(Cl(1,0,0).e1 ≈ Cl(1,1,1).e1)
        alg = Cl(3,0,0)
        @test alg.e1 ≈ alg.e1
        @test !(alg.e1 ≈ alg.e1 + 1e-6)
        @test alg.e1 ≈ alg.e1 + 1e-15
        @test !isapprox(alg.e1, alg.e1 + 1e-6, atol=1e-8)
        @test alg.e1 ≈ (alg.e1 + 1e-6) atol = 2e-6
        @test !isapprox(alg.e1, alg.e1 + 1e-6, rtol=1e-8)
        @test alg.e1 ≈ (alg.e1 + 1e-6) rtol = 2e-6

        @test !isapprox(alg.e1, 1)
        @test !isapprox(1, alg.e2)
        @test isapprox(alg.e1, 1, atol=2)
        @test !isapprox(alg.e1, 1, atol=1)
        @test isapprox(1, alg.e1, atol=2)
        @test !isapprox(1, alg.e1, atol=1)
        @test isapprox(0alg.e1, 0alg.e1)
        @test isapprox(0alg.e1, 0alg.e2)
    end
    @testset "Algebra" begin
        @test CliffordAlgebra(1,0,0) == CliffordAlgebra(1)
        @test CliffordAlgebra(0,1,0) == CliffordAlgebra(0,1)
        @test CliffordAlgebra(0,0,1) == CliffordAlgebra(0,0,1)

        @test_throws AssertionError CliffordAlgebra(-1)
        @test_throws AssertionError CliffordAlgebra(0,-1)
        @test_throws AssertionError CliffordAlgebra(0,0,-1)

        @test_throws AssertionError CliffordAlgebra(0,0,0, (:e1,))
        @test_throws AssertionError CliffordAlgebra(1,0,0, (:e1,:e2,))

        @test typeof(CliffordAlgebra(1,0,0, (:e1,))) <: CliffordAlgebra
        @test typeof(CliffordAlgebra(1,1,0, (:e₊,:e₋))) <: CliffordAlgebra
        @test typeof(CliffordAlgebra(1,1,1, (:e₊,:e₋,:e₀))) <: CliffordAlgebra
        
        @test CliffordAlgebra(:Hyper) == CliffordAlgebra(:Hyperbolic)
        @test CliffordAlgebra(:Complex) == CliffordAlgebra(:ℂ)
        @test CliffordAlgebra(:Dual) == CliffordAlgebra(:Grassmann)
        @test CliffordAlgebra(:Grassmann2D) == CliffordAlgebra(:G2)
        @test CliffordAlgebra(:Grassmann3D) == CliffordAlgebra(:G3)
        @test CliffordAlgebra(:Quaternions) == CliffordAlgebra(:ℍ)
        @test CliffordAlgebra(:Cl2) == CliffordAlgebra(2)
        @test CliffordAlgebra(:Cl3) == CliffordAlgebra(3)
        @test CliffordAlgebra(:Spacetime) == CliffordAlgebra(:STA)
        @test CliffordAlgebra(:PGA2D) == CliffordAlgebra(:Projective2D) == CliffordAlgebra(:Plane2D)
        @test CliffordAlgebra(:PGA3D) == CliffordAlgebra(:Projective3D) == CliffordAlgebra(:Plane3D)
        @test CliffordAlgebra(:CGA2D) == CliffordAlgebra(:Conformal2D)
        @test CliffordAlgebra(:CGA3D) == CliffordAlgebra(:Conformal3D)
        @test CliffordAlgebra(:DCGA3D) == CliffordAlgebra(:DoubleConformal3D)
        @test CliffordAlgebra(:TCGA3D) == CliffordAlgebra(:TripleConformal3D)
        @test CliffordAlgebra(:DCGSTA) == CliffordAlgebra(:DoubleConformalSpacetime)
        @test CliffordAlgebra(:QCGA) == CliffordAlgebra(:QuadricConformal)

        @test_throws ArgumentError CliffordAlgebra(:UnknownAlgebra)

        cl2 = CliffordAlgebra(2)

        @test basesymbol(cl2, 1) == :𝟏
        @test basesymbol(cl2, 2) == :e1
        @test basesymbol(cl2, 3) == :e2
        @test basesymbol(cl2, 4) == :e1e2

        @test character(cl2) == -1
        @test dimension(cl2) == 4
        @test order(cl2) == 2
        @test signature(cl2) == (2,0,0)

        io = IOBuffer()
        signaturetable(io, cl2)
        @test String(take!(io)) == "┌────┬────┐\n│ e1 │ +1 │\n│ e2 │ +1 │\n└────┴────┘\n"
        cayleytable(io, cl2)
        @test String(take!(io)) == "┌───────┬──────────────┬───────┐\n│  +1   │  +e1    +e2  │ +e1e2 │\n├───────┼──────────────┼───────┤\n│  +e1  │  +1    +e1e2 │  +e2  │\n│  +e2  │ -e1e2   +1   │  -e1  │\n├───────┼──────────────┼───────┤\n│ +e1e2 │  -e2    +e1  │  -1   │\n└───────┴──────────────┴───────┘\n"
    end

    @testset "MultiVector" begin
        cl3 = CliffordAlgebra(:Cl3)

        @test scalar(MultiVector(cl3,1)) == 1
        @test MultiVector(cl3, (1,2,3)).e1 == 1
        @test MultiVector(cl3, (1,2,3)).e2 == 2
        @test MultiVector(cl3, (1,2,3)).e3 == 3
        @test MultiVector(cl3, (1,2,3)).𝟏 == 0
        @test MultiVector(cl3, (1,2,3)).e1e2 == 0
        @test MultiVector(cl3, (1,2,3)).e2e3 == 0
        @test MultiVector(cl3, (1,2,3)).e3e1 == 0
        @test MultiVector(cl3, (1,2,3)).e1e2e3 == 0

        @test coefficient(MultiVector(cl3, (1,2,3)), 1) == 0
        @test coefficient(MultiVector(cl3, (1,2,3)), 2) == 1
        @test coefficient(MultiVector(cl3, (1,2,3)), 3) == 2
        @test coefficient(MultiVector(cl3, (1,2,3)), 4) == 3
        @test coefficient(MultiVector(cl3, (1,2,3)), 5) == 0

        @test coefficient(MultiVector(cl3, (1,2,3)), :e1) == 1
        @test coefficient(MultiVector(cl3, (1,2,3)), :e2) == 2
        @test coefficient(MultiVector(cl3, (1,2,3)), :e3) == 3

        @test_throws AssertionError MultiVector(cl3, (1,2))

        mv = 1 + cl3.e1 - cl3.e2 + cl3.e3 - cl3.e1e2 + cl3.e2e3 - cl3.e3e1 + cl3.e1e2e3

        @test MultiVector(cl3, 0) == zero(mv)
        @test MultiVector(cl3, 1) == one(mv)

        @test iszero(zero(mv))
        @test isone(one(mv))

        @test algebra(mv) == cl3

        @test eltype(mv) == typeof(scalar(mv))

        @test convert(Float64, one(mv)) == one(Float64)
        @test_throws InexactError convert(Float64, mv)

        mv2 = cl3.e1 + 0 * cl3.e2

        @test prune(mv) == mv
        @test typeof(prune(mv2)) != typeof(mv2)

        @test isgrade(mv2, 1)
        @test !isgrade(mv2, 0)
        @test !isgrade(mv2, 2)
        @test !isgrade(mv2, 3)

        @test !isgrade(mv, 0)
        @test !isgrade(mv, 1)
        @test !isgrade(mv, 2)
        @test !isgrade(mv, 3)

        @test isgrade(grade(mv,2),2)

        @test even(mv2) == zero(mv2)
        @test odd(mv2) == mv2

        @test grin(mv2) == -mv2

        @test dual(mv2) == cl3.e2e3

        @test reverse(cl3.𝟏) == cl3.𝟏
        @test reverse(cl3.e2e3) == -cl3.e2e3
        @test reverse(cl3.e1e2e3) == -cl3.e1e2e3

        @test reverse(mv) == ~mv

        @test pseudoscalar(cl3) == cl3.e1e2e3

        mva = cl3.e1 - cl3.e1e2
        mvb = 1 + cl3.e2 + cl3.e1e2e3

        @assert matrix(mva) * matrix(mvb) == matrix(mva*mvb)
        @assert matrix(mva) * vector(mvb) == vector(mva*mvb)
    end

    @testset "arithmetic for $cl" for cl in (
            CliffordAlgebra(3,0,0),
            CliffordAlgebra(2,1,0),
            CliffordAlgebra(1,2,0),
            CliffordAlgebra(2,0,1),
            CliffordAlgebra(1,0,2),
            CliffordAlgebra(0,2,1),
            CliffordAlgebra(0,1,2),
            CliffordAlgebra(1,1,1)
        ) 

        e1 = cl.e1
        e2 = cl.e2
        e3 = cl.e3
        e12 = cl.e1e2
        e23 = cl.e2e3
        e31 = cl.e3e1
        e123 = cl.e1e2e3

        s1 = scalar(e1*e1)
        s2 = scalar(e2*e2)
        s3 = scalar(e3*e3)

        @test character(cl) == scalar(e123*e123)

        @testset "inference" begin
            @inferred 3*e1
            @inferred 3.0*e1
            @inferred e1*3f0
            @inferred (3e1) * (e2)
            @inferred e1 ∧ 2.3
            @inferred -e1
        
            @inferred dual(e1)
            @inferred dual(e1 + e12)
            @inferred e1 ∨ (e23 + 1)
        
            @inferred scalar(e1)
            @inferred scalar(3*e1)
            @inferred reverse(e1 + 3*e12)
            @inferred norm(e1 + 3*e12)
        end

        @testset "addition/subtraction" begin
            @test e1 + e2 == e2 + e1
            @test e1 - e2 == -(e2 - e1)
            @test e3 + 1 == 1 + e3
        end

        @testset "geometric product" begin 
            @test e1 * e1 == s1
            @test e2 * e2 == s2
            @test e3 * e3 == s3

            @test e1 * e2 == e12
            @test e2 * e3 == e23
            @test e1 * e3 == -e31
            @test e1 * e2 * e3 == e123
        end

        @testset "exterior product" begin 
            @test e1 ∧ e1 == 0
            @test e2 ∧ e2 == 0
            @test e3 ∧ e3 == 0

            @test e1 ∧ e2 == e12
            @test e2 ∧ e3 == e23
            @test e3 ∧ e1 == e31

            @test e12 ∧ e31 == 0
            @test e23 ∧ e31 == 0
            @test e12 ∧ e23 == 0
            
            @test e1 ∧ e23 == e123
            @test e12 ∧ e3 == e123

            @test 2 ∧ e1 == MultiVector(cl,2) ∧ e1
            @test 2 ∧ 3 == MultiVector(cl,2) ∧ MultiVector(cl,3)
        end

        @testset "fatdot product" begin 
            @test e1 ⋅ e1 == s1
            @test e2 ⋅ e2 == s2
            @test e3 ⋅ e3 == s3

            @test e1 ⋅ e12 == s1*e2
            @test e1 ⋅ e23 == 0
            @test e1 ⋅ e31 == -s1*e3

            @test e123 ⋅ e1 == s1*e23
            @test e123 ⋅ e2 == s2*e31
            @test e123 ⋅ e3 == s3*e12
            @test e123 ⋅ e12 == -s1*s2*e3
            @test e123 ⋅ e23 == -s2*s3*e1
            @test e123 ⋅ e31 == -s3*s1*e2
            @test e123 ⋅ e123 == -s1*s2*s3

            @test 2 ⋅ e1 == MultiVector(cl,2) ⋅ e1
            @test 2 ⋅ 3 == MultiVector(cl,2) ⋅ MultiVector(cl,3) 
        end

        @testset "left contraction product" begin
            @test e1 ⨼ e1 == s1
            @test e2 ⨼ e2 == s2
            @test e3 ⨼ e3 == s3

            @test e1 ⨼ e12 == s1*e2
            @test e1 ⨼ e23 == 0
            @test e1 ⨼ e31 == -s1*e3

            @test e123 ⨼ e1 == 0
            @test e123 ⨼ e2 == 0
            @test e123 ⨼ e3 == 0

            @test e123 ⨼ e12 == 0
            @test e123 ⨼ e23 == 0
            @test e123 ⨼ e31 == 0

            @test e123 ⨼ e123 == -s1*s2*s3
            
            @test 2 ⨼ e1 == MultiVector(cl,2) ⨼ e1
            @test 2 ⨼ 3 == MultiVector(cl,2) ⨼ MultiVector(cl,3)
        end

        @testset "right contraction product" begin
            @test e1 ⨽ e1 == s1
            @test e2 ⨽ e2 == s2
            @test e3 ⨽ e3 == s3

            @test e1 ⨽ e12 == 0
            @test e1 ⨽ e23 == 0
            @test e1 ⨽ e31 == 0

            @test e123 ⨽ e1 == s1*e23
            @test e123 ⨽ e2 == s2*e31
            @test e123 ⨽ e3 == s3*e12

            @test e123 ⨽ e12 == -s1*s2*e3
            @test e123 ⨽ e23 == -s2*s3*e1
            @test e123 ⨽ e31 == -s3*s1*e2

            @test e123 ⨽ e123 == -s1*s2*s3
            
            @test 2 ⨽ e1 == MultiVector(cl,2) ⨽ e1
            @test 2 ⨽ 3 == MultiVector(cl,2) ⨽ MultiVector(cl,3)
        end

        @testset "scalar product" begin
            @test e1 ⋆ e1 == s1
            @test e2 ⋆ e2 == s2
            @test e3 ⋆ e3 == s3

            @test e1 ⋆ e12 == 0
            @test e1 ⋆ e23 == 0
            @test e1 ⋆ e31 == 0

            @test e123 ⋆ e1 == 0
            @test e123 ⋆ e2 == 0
            @test e123 ⋆ e3 == 0

            @test e123 ⋆ e123 == -s1 * s2 * s3

            @test 2 ⋆ e1 == MultiVector(cl,2) ⋆ e1
            @test 2 ⋆ 3 == MultiVector(cl,2) ⋆ MultiVector(cl,3)
        end

        @testset "regressive product" begin
            @test e1 ∨ e1 == 0
            @test e2 ∨ e2 == 0
            @test e3 ∨ e3 == 0

            @test e1 ∨ e23 == 1
            @test e2 ∨ e31 == 1
            @test e3 ∨ e12 == 1

            @test e123 ∨ e1 == e1
            @test e123 ∨ e2 == e2
            @test e123 ∨ e3 == e3

            @test e123 ∨ e12 == e12
            @test e123 ∨ e23 == e23
            @test e123 ∨ e31 == e31

            @test e123 ∨ e123 == e123

            @test 2 ∨ e1 == MultiVector(cl,2) ∨ e1
            @test 2 ∨ 3 == MultiVector(cl,2) ∨ MultiVector(cl,3)

            @test 1 ∨ 2 == 0
            @test e1 ∨ 2 == 0
            @test e12 ∨ 2 == 0
            @test e123 ∨ 2 == 2
        end

        @testset "commutator product" begin
            @test 1 ×₋ 2 == 0
            @test 2 ×₋ e1 == 0 
            @test 2 ×₋ e2 == 0 
            @test 2 ×₋ e3 == 0 
            @test 2 ×₋ e12 == 0 
            @test 2 ×₋ e23 == 0
            @test 2 ×₋ e31 == 0 
            @test 2 ×₋ e123 == 0 

            @test e1 ×₋ e1 == 0
            @test e2 ×₋ e2 == 0
            @test e3 ×₋ e3 == 0

            @test e1 ×₋ e2 == e12
            @test e1 ×₋ e3 == -e31
            @test e2 ×₋ e3 == e23

            @test e12 ×₋ e12 == 0
            @test e23 ×₋ e23 == 0
            @test e31 ×₋ e31 == 0

            @test e23 ×₋ e12 == s2 * e31
            @test e23 ×₋ e31 == -s3 * e12
            @test e12 ×₋ e31 == s1 * e23

            @test e123 ×₋ e123 == 0
        end

        @testset "anti-commutator product" begin
            @test 1 ×₊ 2 == 2
            @test 2 ×₊ e1 == 2*e1 
            @test 2 ×₊ e2 == 2*e2 
            @test 2 ×₊ e3 == 2*e3
            @test 2 ×₊ e12 == 2*e12 
            @test 2 ×₊ e23 == 2*e23
            @test 2 ×₊ e31 == 2*e31 
            @test 2 ×₊ e123 == 2*e123 

            @test e1 ×₊ e1 == s1
            @test e2 ×₊ e2 == s2
            @test e3 ×₊ e3 == s3

            @test e1 ×₊ e2 == 0
            @test e1 ×₊ e3 == 0
            @test e2 ×₊ e3 == 0

            @test e12 ×₊ e12 == -s1 * s2
            @test e23 ×₊ e23 == -s2 * s3
            @test e31 ×₊ e31 == -s3 * s1

            @test e23 ×₊ e12 == 0
            @test e23 ×₊ e31 == 0
            @test e12 ×₊ e31 == 0

            @test e123 ×₊ e123 == -s1 * s2 * s3
        end

        @testset "sandwich product" for n = 1:100
            r = -100:100
            mva = rand(r) + rand(r) * e1 + 
                  rand(r) * e2 + rand(r) * e3 + 
                  rand(r) * e12 + rand(r) * e23 +
                  rand(r) * e31 + rand(r) * e123

            mvb = rand(r) + rand(r) * e1 + 
                  rand(r) * e2 + rand(r) * e3 + 
                  rand(r) * e12 + rand(r) * e23 +
                  rand(r) * e31 + rand(r) * e123

            @test mva ≀ mvb == mva * mvb * reverse(mva)
        end

        @testset "coefficients" begin
            pga = typeof(CliffordAlgebra(:PGA3D))
            mv_int = MultiVector(pga, (1:16...,), (101:116...,))
            mv_float = MultiVector(pga, (1:16...,), (101.5:116.5...,))
            
            fib_ntuple = (1, 2, 3, 5, 8, 13)
            fib_vec = [1, 2, 3, 5, 8, 13]
            sym_ntuple = (:𝟏, :e1, :e2, :e0, :e2e3, :e1e2e0)
            sym_vec = [:𝟏, :e1, :e2, :e0, :e2e3, :e1e2e0]
            
            @test coefficients(mv_int, fib_ntuple) == fib_ntuple .+ 100
            @test coefficients(mv_int, fib_vec) == fib_vec .+ 100
            @test coefficients(mv_float, fib_ntuple) == fib_ntuple .+ 100.5
            @test coefficients(mv_float, fib_vec) == fib_vec .+ 100.5
            
            @test coefficients(mv_int, sym_ntuple) == fib_ntuple .+ 100
            @test coefficients(mv_int, sym_vec) == fib_vec .+ 100
            @test coefficients(mv_float, sym_ntuple) == fib_ntuple .+ 100.5
            @test coefficients(mv_float, sym_vec) == fib_vec .+ 100.5
        end
        
        @testset "broadcasted .*" begin
            pga = typeof(CliffordAlgebra(:PGA2D))
            mvs = [
                MultiVector(pga, (1:3...,), (2,2,2)),
                MultiVector(pga, (1:6...,), (3,3,3,3,3,3)),
                MultiVector(pga, (4:8...,), (5,5,5,5,5)),
                MultiVector(pga, (1:8...,), (7,7,7,7,7,7,7,7)),
            ]
            
            for mv1 in mvs, mv2 in mvs
                @test isapprox(vector(mv1 .* mv2), vector(mv1) .* vector(mv2))
            end
        end

        @testset "broadcasted ./" begin
            pga = typeof(CliffordAlgebra(:PGA2D))

            function semi_safe_divide(x, y)
                if x == 0 && y == 0
                    return 0
                else
                    return x / y
                end
            end

            mv_zero = MultiVector(pga, (1:3...,), (0, 0, 0))
            mv_small = MultiVector(pga, (1:3...,), (2, 2, 2))
            mv_full = MultiVector(pga, (1:8...,), (7, 7, 7, 7, 7, 7, 7, 7))

            @test isapprox(vector(mv_small ./ mv_zero), [Inf, Inf, Inf, 0, 0, 0, 0, 0])
            @test isapprox(vector(mv_small ./ mv_zero), semi_safe_divide.(vector(mv_small), vector(mv_zero)))
            @test_throws AssertionError mv_full ./ mv_zero

            @test isapprox(vector(mv_zero ./ mv_small), semi_safe_divide.(vector(mv_zero), vector(mv_small)))
            @test isapprox(vector(mv_small ./ mv_small), semi_safe_divide.(vector(mv_small), vector(mv_small)))
            @test_throws AssertionError mv_full ./ mv_small

            @test isapprox(vector(mv_zero ./ mv_full), semi_safe_divide.(vector(mv_zero), vector(mv_full)))
            @test isapprox(vector(mv_small ./ mv_full), semi_safe_divide.(vector(mv_small), vector(mv_full)))
            @test isapprox(vector(mv_full ./ mv_full), semi_safe_divide.(vector(mv_full), vector(mv_full)))
        end

        @testset "norm" begin
            r = -100:100
            mv  = rand(r) + rand(r) * e1 + 
                  rand(r) * e2 + rand(r) * e3 + 
                  rand(r) * e12 + rand(r) * e23 +
                  rand(r) * e31 + rand(r) * e123
            
            @test norm(mv) ≈ sqrt(Complex(scalar(mv*reverse(mv))))
        end

        @testset "misc functions" begin
            mv = 1 + e1 - e2 + e3 + e12 - e23 + e31 - e123
            @test polarize(mv)*pseudoscalar(algebra(mv)) == character(algebra(mv)) * mv
            @test Λᵏ(mv,0) == grade(mv,0)
            @test Λᵏ(mv,1) == grade(mv,1)
            @test Λᵏ(mv,2) == grade(mv,2)
            @test Λᵏ(mv,3) == grade(mv,3)

            @test maxgrade(mv) == (grade(mv,3),3)
            @test mingrade(mv) == (grade(mv,0),0)

            @test norm_sqr(e1) == s1
            @test norm_sqr(2*e2) == 4*s2
            @test norm_sqr(3*e3) == 9*s3
            @test norm_sqr(e12) == s1*s2
            @test norm_sqr(e23) == s2*s3
            @test norm_sqr(e31) == s3*s1
            @test norm_sqr(e123) == s1*s2*s3

            @test norm(2*e1)  == 2*sqrt(Complex(s1))
            @test norm(-2*e1) == 2*sqrt(Complex(s1))

            @test norm(e1+2*e31) ≈ sqrt(norm(e1)^2 + 4*norm(e31)^2)

            if !iszero(s1)
                @test inv(e1) * e1 == e1 * inv(e1) == 1 
            else
                @test_throws SingularException inv(e1)
                @test_throws SingularException inv(e12)
                @test_throws SingularException inv(e31)
                @test_throws SingularException inv(e123)
            end

            if !iszero(s2)
                @test inv(e2) * e2 == e2 * inv(e2) == 1
            else
                @test_throws SingularException inv(e2)
                @test_throws SingularException inv(e12)
                @test_throws SingularException inv(e23)
                @test_throws SingularException inv(e123)
            end

            if !iszero(s3)
                @test inv(e3) * e3 == e3 * inv(e3) == 1
            else
                @test_throws SingularException inv(e3)
                @test_throws SingularException inv(e31)
                @test_throws SingularException inv(e23)
                @test_throws SingularException inv(e123)
            end

            if !iszero(s1) && !iszero(s2)
                @test inv(e12) * e12 == e12 * inv(e12) == 1
            end

            if !iszero(s1) && !iszero(s3)
                @test inv(e31) * e31 == e31 * inv(e31) == 1
            end

            if !iszero(s2) && !iszero(s3)
                @test inv(e23) * e23 == e23 * inv(e23) == 1
            end

            if !iszero(s1) && !iszero(s2) && !iszero(s3)
                @test (e1+e2+e3) \ (e1+e2+e3) ≈ 1
                @test (e12 + e23 + e31) / (e12 + e23 + e31) ≈ 1
                @test inv(e123) * e123 == 1
            end

            for mv in (e1, e2, e3, e12, e23, e31, e123)
                @inferred exp(mv)
                if scalar(mv*mv) < 0
                    @test exp(mv) == cos(1) + sin(1) * mv
                elseif scalar(mv*mv) > 0
                    @test exp(mv) == cosh(1) + sinh(1) * mv
                else
                    @test exp(mv) == 1 + mv
                end
                @test convert(Float32,exp(mv) * exp(-mv)) ≈ 1
                @test exp(0*mv) == 1
            end

            @test exp(1 + e1 + e12 + e123) ≈ exp(1 + e123) * exp(e1 + e12)

            mv = 1 + e1 - e2 - e3 + e12 - e23 + e31 - e123
            
            M = [1 0 0; 0 1 0 ; 0 0 1]
            @test outermorphism(M, mv) == mv

            M = [0 0 0; 0 0 0 ; 0 0 0]
            @test outermorphism(M, mv) == 1

            M = [0 1 0; 1 0 0 ; 0 0 1]
            @test outermorphism(M, mv) == 1 + e2 - e1 - e3 - e12 + e31 - e23 + e123

            M = [0 0 1; 0 1 0 ; 1 0 0]
            @test outermorphism(M, mv) == 1 + e3 - e2 - e1 - e23 + e12 - e31 + e123

            M = [0 0 1; 0 1 0 ; 0 0 0]
            @test outermorphism(M, mv) == 1 - e2 - e1 + e12

            @test !(e1 ≈ e2)
            @test e1 ≈ e1
            @test !(e1 ≈ 1.00001*e1)
            @test e1 ≈ 1.00001*e1 atol = 0.001
            @test !(e2 ≈ e2 + 1e-6*e12)
            @test (e2 ≈ e2 + 1e-8*e12)
        end

    end
end
