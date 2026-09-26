using Test
using Random
import AbstractAlgebra as AA
using QuantumEnvelopingAlgebras
using QuantumEnvelopingAlgebras: UqSl2Elem, tensor_terms, QuantumPlaneElem, PlaneMonomial
include("reference_word_rewriter.jl")
include("test_representations.jl")
include("test_hopf.jl")
include("test_specialization.jl")
include("test_quantum_plane.jl")

function check_algebra(R, q)
    U, (E, F, K, Ki) = uqsl2(R, q)
    d = q - inv(q)
    @test AA.parent(E) === U
    @test AA.base_ring(U) === R
    @test AA.ngens(U) == 4
    @test quantum_parameter(U) == q
    @test K*Ki == Ki*K == one(U)
    @test K*E*Ki == q^2*E
    @test K*F*Ki == q^(-2)*F
    @test commutator(E, F) == (K - Ki)/d
    @test E*F == F*E + (K-Ki)/d
    @test K^(-3) == Ki^3
    @test E^0 == one(U)
    @test zero(U)^0 == one(U)
    @test inv(2*K^3)*2*K^3 == one(U)
    @test K^Int8(-3) == Ki^3
    @test_throws OverflowError K^typemin(Int8)
    @test_throws UnsupportedInverseError inv(E)
    @test_throws UnsupportedInverseError inv(zero(U))
    @test pbw_coefficient(E*F, 1, 0, 1) == one(R)
    @test pbw_coefficient(E*F, 0, 1, 0) == inv(d)
    @test pbw_coefficient(E*F, 0, -1, 0) == -inv(d)
    @test iszero(commutator(casimir(U), E))
    @test iszero(commutator(casimir(U), F))
    @test iszero(commutator(casimir(U), K))

    x = E*F + K^(-2) + q*F
    @test x == normal_form(x)
    @test hash(x) == hash(copy(x)) == hash(deepcopy(x))
    @test sprint(show, x) == sprint(show, copy(x))
    @test AA.parent(deepcopy(x)) === U
    @test x == copy(x)
    @test isequal(x, copy(x))
    @test x == U(x)
    @test_throws ArgumentError x.terms
    raw = Dict(PBWMonomial(0, 0, 1) => one(R))
    from_raw = UqSl2Elem(U, raw)
    empty!(raw)
    @test from_raw == E
    @test U(2) == 2
    @test !isequal(U(2), 2)
    @test_throws DomainError basis_monomial(U, -1, 0, 0)
    @test_throws OverflowError basis_monomial(U, big(typemax(Int))+1, 0, 0)
    @test_throws DomainError qinteger(q, -1)
    @test_throws OverflowError pbw_terms(basis_monomial(U, typemax(Int), 0, typemax(Int)))

    V, (e, f, k, ki) = uqsl2(R, q)
    @test E != e
    @test_throws ParentMismatchError E + e
    @test_throws ParentMismatchError E * e
    @test_throws ParentMismatchError U(e)

    rng = MersenneTwister(1229)
    atoms = [E, F, K, Ki, one(U), -one(U)]
    for _ in 1:30
        x, y, z = (sum(rand(rng, atoms) for _ in 1:3) for _ in 1:3)
        @test (x*y)*z == x*(y*z)
        @test x*(y+z) == x*y + x*z
        @test (x+y)*z == x*z + y*z
        @test iszero(x - x)
    end
end

@testset "U_q(sl_2)" begin
    R, q = AA.rational_function_field(AA.QQ, :q)
    check_algebra(R, q)
    check_reference_rewriter(uqsl2(R, q)[1], q)
    check_representations(uqsl2(R, q)[1], q)
    check_hopf(R, q)
    check_algebra(AA.QQ, AA.QQ(2))
    check_reference_rewriter(uqsl2(AA.QQ, AA.QQ(2))[1], AA.QQ(2))
    check_representations(uqsl2(AA.QQ, AA.QQ(2))[1], AA.QQ(2))
    check_hopf(AA.QQ, AA.QQ(2))
    q0 = AA.QQ(3)/AA.QQ(2)
    check_algebra(AA.QQ, q0)
    check_reference_rewriter(uqsl2(AA.QQ, q0)[1], q0)
    check_representations(uqsl2(AA.QQ, q0)[1], q0)
    check_hopf(AA.QQ, q0)
    @test_throws DomainError uqsl2(AA.QQ, 0)
    @test_throws DomainError uqsl2(AA.QQ, 1)
    @test_throws DomainError uqsl2(AA.QQ, -1)
    @test_throws ArgumentError uqsl2(AA.ZZ, 2)
    @test_throws ArgumentError uqsl2(AA.QQ, 2; max_terms=0)
    limited, (e, f, k, ki) = uqsl2(AA.QQ, 2; max_terms=2)
    @test e + f == f + e
    @test_throws ComputationLimitError e*f
    @test_throws ComputationLimitError e + f + k
    check_specialization()
end
