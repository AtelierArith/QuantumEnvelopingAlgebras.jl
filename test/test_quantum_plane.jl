function check_quantum_plane(R, q)
    P, (X,Y) = quantum_plane(R,q)
    @test AA.parent(X) === P
    @test AA.base_ring(P) === R
    @test AA.ngens(P) == 2
    @test quantum_parameter(P) == q
    @test X*Y == q*Y*X
    @test X^3*Y^4 == q^12*Y^4*X^3
    @test basis_monomial(P,2,3) == Y^2*X^3
    @test pbw_coefficient(X*Y,1,1) == q
    @test pbw_coefficient(X*Y,0,0) == zero(R)
    @test P(2) == 2
    @test !isequal(P(2),2)
    @test inv(2*one(P)) == P(inv(R(2)))
    @test_throws UnsupportedInverseError inv(X)
    @test_throws UnsupportedInverseError inv(zero(P))
    @test zero(P)^0 == one(P)
    @test X^0 == one(P)
    @test P(2)^(-2) == P(inv(R(4)))
    @test_throws OverflowError X^typemin(Int8)

    z = (X+Y)^4 + q*Y^2*X + one(P)
    @test z == normal_form(z)
    @test hash(z) == hash(copy(z)) == hash(deepcopy(z))
    @test sprint(show,z) == sprint(show,deepcopy(z))
    @test AA.parent(deepcopy(z)) === P
    @test_throws ArgumentError z.terms
    raw = Dict(PlaneMonomial(1,0) => one(R))
    snapshot = QuantumPlaneElem(P,raw)
    empty!(raw)
    @test snapshot == Y
    @test_throws DomainError basis_monomial(P,-1,0)
    @test_throws DomainError basis_monomial(P,0,-1)
    @test_throws OverflowError basis_monomial(P,big(typemax(Int))+1,0)
    @test_throws OverflowError pbw_terms(basis_monomial(P,typemax(Int),typemax(Int)))
    @test_throws OverflowError basis_monomial(P,typemax(Int),0)*Y

    Q, (x,y) = quantum_plane(R,q)
    @test X != x
    @test_throws ParentMismatchError X + x
    @test_throws ParentMismatchError X * x
    @test_throws ParentMismatchError P(x)

    # Count inversions in a free word, independently of the package product.
    for n in 0:7, letters in Iterators.product(ntuple(_ -> (:X,:Y),n)...)
        word = Tuple(letters)
        inversions = count(i < j && word[i] == :X && word[j] == :Y
            for i in 1:length(word) for j in i+1:length(word))
        expected = q^inversions*basis_monomial(P,count(==(:Y),word),count(==(:X),word))
        actual = foldl(*,(letter == :X ? X : Y for letter in word); init=one(P))
        @test actual == expected
    end

    rng = MersenneTwister(733)
    atoms = (X,Y,one(P),-one(P))
    for _ in 1:30
        a,b,c = (sum(rand(rng,atoms) for _ in 1:3) for _ in 1:3)
        @test (a*b)*c == a*(b*c)
        @test a*(b+c) == a*b+a*c
        @test (a+b)*c == a*c+b*c
    end
end

@testset "Quantum plane" begin
    R,q = AA.rational_function_field(AA.QQ,:q)
    check_quantum_plane(R,q)
    check_quantum_plane(AA.QQ,AA.QQ(2))
    check_quantum_plane(AA.QQ,AA.QQ(1))
    @test_throws DomainError quantum_plane(AA.QQ,0)
    @test_throws ArgumentError quantum_plane(AA.ZZ,2)
    @test_throws ArgumentError quantum_plane(AA.QQ,2; max_terms=0)
    limited,(X,Y) = quantum_plane(AA.QQ,2; max_terms=2)
    @test X+Y == Y+X
    @test_throws ComputationLimitError (X+Y)^2
end
