function check_specialization()
    R, q = AA.rational_function_field(AA.QQ, :q)
    U, (E,F,K,Ki) = uqsl2(R, q)
    for point in (AA.QQ(2), AA.QQ(3)/AA.QQ(2))
        phi = specialization_map(U, AA.QQ, point)
        @test AA.parent(phi(E)) === phi.target
        @test quantum_parameter(phi.target) == point
        @test phi(E*F) == phi(E)*phi(F)
        @test phi(E+F+K+Ki) == phi(E)+phi(F)+phi(K)+phi(Ki)
        a = ((q+1)/(q-3))*E + inv(q)*Ki
        b = ((q-4)/(q+3))*F + K
        @test phi(a*b) == phi(a)*phi(b)
        @test phi(a+b) == phi(a)+phi(b)
        @test phi(casimir(U)) == casimir(phi.target)
        @test phi(U((q-point)/(q-point))) == one(phi.target)
        @test_throws SpecializationError phi(U(inv(q-point)))
        @test_throws ParentMismatchError phi(AA.gens(uqsl2(R,q)[1])[1])
    end
    @test_throws DomainError specialization_map(U, AA.QQ, 0)
    @test_throws DomainError specialization_map(U, AA.QQ, 1)
    @test_throws DomainError specialization_map(U, AA.QQ, -1)
    @test_throws ArgumentError specialization_map(uqsl2(AA.QQ,2)[1], AA.QQ, 2)

    U2, (e,f,k,ki) = uqsl2(R, q^2)
    phi2 = specialization_map(U2, AA.QQ, 2)
    @test quantum_parameter(phi2.target) == AA.QQ(4)
    @test phi2(e*f-f*e) == (phi2(k)-phi2(ki))/(AA.QQ(4)-inv(AA.QQ(4)))

    bad, _ = uqsl2(R, inv(q-2))
    @test_throws SpecializationError specialization_map(bad, AA.QQ, 2)
end
