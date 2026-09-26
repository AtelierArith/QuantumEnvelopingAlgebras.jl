struct CoproductMap{T<:UqSl2TensorPower{2}}
    source::UqSl2Algebra
    target::T
end

"""
    coproduct_map(U, T2)

Return the algebra map `Delta: U -> T2` (with `T2 == tensor_power(U, 2)`) given
by `K -> K ⊗ K`, `E -> E ⊗ 1 + K ⊗ E`, `F -> F ⊗ K^-1 + 1 ⊗ F`.
"""
function coproduct_map(U::UqSl2Algebra, T::UqSl2TensorPower{2})
    T.algebra === U || throw(ParentMismatchError())
    return CoproductMap(U, T)
end

function (map::CoproductMap)(x::UqSl2Elem)
    parent(x) === map.source || throw(ParentMismatchError())
    U = map.source
    T = map.target
    E, F, K, Ki = gens(U)
    unit = one(U)
    delta_E = tensor(T, E, unit) + tensor(T, K, E)
    delta_F = tensor(T, F, Ki) + tensor(T, unit, F)
    result = zero(T)
    for (m, coeff) in _terms(x)
        delta_Kb = tensor(T, K^m.b, K^m.b)
        result += coeff*(delta_F^m.a * delta_Kb * delta_E^m.c)
    end
    return result
end

"Apply the coproduct to one slot of a tensor element, flattening into T."
function coproduct_on_slot(T::UqSl2TensorPower{M}, x::UqSl2TensorElem{N}, slot::Integer) where {M,N}
    M == N + 1 || throw(ArgumentError("target tensor power must have one more slot"))
    1 <= slot <= N || throw(BoundsError(1:N, slot))
    parent(x).algebra === T.algebra || throw(ParentMismatchError())
    U = T.algebra
    T2 = tensor_power(U, 2)
    delta = coproduct_map(U, T2)
    result = zero(T)
    for (key, coeff) in _tensor_terms(x)
        selected = key[slot]
        image = delta(basis_monomial(U, selected.a, selected.b, selected.c))
        for (pair, image_coeff) in _tensor_terms(image)
            factors = ntuple(i -> begin
                m = i < slot ? key[i] : i == slot ? pair[1] : i == slot+1 ? pair[2] : key[i-1]
                basis_monomial(U, m.a, m.b, m.c)
            end, M)
            result += (coeff*image_coeff)*tensor(T, factors...)
        end
    end
    return result
end

"""
    counit(x)

Return the counit `epsilon(x)` in the base ring (`K -> 1`, `E, F -> 0`).
"""
function counit(x::UqSl2Elem)
    result = zero(base_ring(x))
    for (m, coeff) in _terms(x)
        m.a == 0 && m.c == 0 && (result += coeff)
    end
    return result
end

"""
    antipode(x)

Return the antipode `S(x)`, an antialgebra map (`K -> K^-1`, `E -> -K^-1 E`,
`F -> -F K`).
"""
function antipode(x::UqSl2Elem)
    U = parent(x)
    E, F, K, Ki = gens(U)
    result = zero(U)
    se = -Ki*E
    sf = -F*K
    for (m, coeff) in _terms(x)
        result += coeff*(se^m.c * K^Base.Checked.checked_neg(m.b) * sf^m.a)
    end
    return result
end
