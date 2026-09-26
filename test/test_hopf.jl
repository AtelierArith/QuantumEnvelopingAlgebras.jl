function tensor_counit(x, slot)
    U = AA.parent(x).algebra
    result = zero(U)
    for (key, coeff) in tensor_terms(x)
        selected = key[slot]
        other = key[3-slot]
        selected_element = basis_monomial(U, selected.a, selected.b, selected.c)
        other_element = basis_monomial(U, other.a, other.b, other.c)
        result += coeff*counit(selected_element)*other_element
    end
    return result
end

function tensor_antipode_product(x, on_left)
    U = AA.parent(x).algebra
    result = zero(U)
    for (key, coeff) in tensor_terms(x)
        left, right = (basis_monomial(U, m.a, m.b, m.c) for m in key)
        result += coeff*(on_left ? antipode(left)*right : left*antipode(right))
    end
    return result
end

function check_hopf(R, q)
    U, (E,F,K,Ki) = uqsl2(R, q)
    T2 = tensor_power(U, 2)
    T3 = tensor_power(U, 3)
    Delta = coproduct_map(U, T2)
    unit = one(U)
    @test AA.parent(Delta(E)) === T2
    @test AA.base_ring(T2) === R
    @test Delta(E) == tensor(T2,E,unit) + tensor(T2,K,E)
    @test Delta(F) == tensor(T2,F,Ki) + tensor(T2,unit,F)
    @test Delta(K) == tensor(T2,K,K)
    @test Delta(Ki) == tensor(T2,Ki,Ki)
    @test Delta(K)*Delta(Ki) == one(T2)
    @test tensor(T2,E,unit)*tensor(T2,unit,F) == tensor(T2,unit,F)*tensor(T2,E,unit)
    @test tensor(T2,E,unit)*tensor(T2,F,unit) == tensor(T2,E*F,unit)
    @test Delta(K)*Delta(E)*Delta(Ki) == q^2*Delta(E)
    @test Delta(K)*Delta(F)*Delta(Ki) == q^(-2)*Delta(F)
    @test Delta(E)*Delta(F) - Delta(F)*Delta(E) == (Delta(K)-Delta(Ki))/(q-inv(q))
    @test counit(E) == counit(F) == zero(R)
    @test counit(K) == counit(Ki) == one(R)
    @test antipode(E) == -Ki*E
    @test antipode(F) == -F*K
    @test antipode(K) == Ki
    @test antipode(Ki) == K

    x = E*F + Ki*E + F^2*K^(-2) + q
    y = F*E + K*F + E^2
    @test Delta(x*y) == Delta(x)*Delta(y)
    @test Delta(x+y) == Delta(x)+Delta(y)
    @test counit(x*y) == counit(x)*counit(y)
    @test antipode(x*y) == antipode(y)*antipode(x)
    @test coproduct_on_slot(T3, Delta(x), 1) == coproduct_on_slot(T3, Delta(x), 2)
    @test tensor_counit(Delta(x), 1) == x
    @test tensor_counit(Delta(x), 2) == x
    @test tensor_antipode_product(Delta(x), true) == counit(x)*unit
    @test tensor_antipode_product(Delta(x), false) == counit(x)*unit
    for generator in (E, F, K, Ki)
        @test coproduct_on_slot(T3, Delta(generator), 1) ==
            coproduct_on_slot(T3, Delta(generator), 2)
        @test tensor_antipode_product(Delta(generator), true) == counit(generator)*unit
        @test tensor_antipode_product(Delta(generator), false) == counit(generator)*unit
    end
    rng = MersenneTwister(2026)
    atoms = (E, F, K, Ki, unit)
    for _ in 1:8
        a = sum(rand(rng, atoms) for _ in 1:3)
        b = sum(rand(rng, atoms) for _ in 1:3)
        @test Delta(a*b) == Delta(a)*Delta(b)
        @test antipode(a*b) == antipode(b)*antipode(a)
    end
    @test hash(Delta(x)) == hash(copy(Delta(x))) == hash(deepcopy(Delta(x)))
    @test AA.parent(deepcopy(Delta(x))) === T2
    @test_throws ArgumentError Delta(x).terms
    @test_throws ParentMismatchError tensor(T2, E, one(uqsl2(R,q)[1]))
    @test_throws ParentMismatchError Delta(E) + tensor(tensor_power(U,2), E, unit)
    @test_throws ParentMismatchError Delta(E) * tensor(T3, E, unit, unit)
    @test_throws ArgumentError tensor_power(U, 0)
    @test_throws BoundsError coproduct_on_slot(T3, Delta(E), 3)
    @test_throws ArgumentError coproduct_on_slot(T2, Delta(E), 1)
end
