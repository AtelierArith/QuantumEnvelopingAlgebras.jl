"A partially defined coefficient evaluation map into a fixed target algebra."
struct SpecializationMap{U,V,C}
    source::U
    target::V
    point::C
end

function _evaluate_polynomial_at(poly, S::AA.Field, point)
    value = zero(S)
    for i in AA.degree(poly):-1:0
        value = value*point + S(AA.coeff(poly, i))
    end
    return value
end

function _evaluate_rational_at(coeff, S::AA.Field, point)
    numerator_value = _evaluate_polynomial_at(numerator(coeff), S, point)
    denominator_value = _evaluate_polynomial_at(denominator(coeff), S, point)
    iszero(denominator_value) && throw(SpecializationError("coefficient denominator vanishes at the specialization point"))
    return numerator_value*inv(denominator_value)
end

"""
    specialization_map(U, S, point)

Return the coefficient evaluation map from `U` (over `QQ(q)`) to the algebra
over the exact characteristic-zero field `S` at `q = point`. A coefficient whose
denominator vanishes at `point` throws [`SpecializationError`](@ref). The result
`phi` exposes its target algebra as `phi.target`.
"""
function specialization_map(U::UqSl2Algebra, S::AA.Field, point)
    R = base_ring(U)
    R isa AA.Generic.RationalFunctionField && base_ring(R) === AA.QQ ||
        throw(ArgumentError("source coefficients must belong to a univariate rational function field over QQ"))
    is_exact_type(elem_type(S)) && iszero(characteristic(S)) ||
        throw(ArgumentError("target must be an exact characteristic-zero field"))
    q0 = try
        S(point)
    catch
        throw(ArgumentError("specialization point cannot be coerced into the target field"))
    end
    iszero(q0) && throw(DomainError(q0, "specialization point must be nonzero"))
    iszero(q0-inv(q0)) && throw(DomainError(q0, "specialization point cannot be 1 or -1"))
    target_q = _evaluate_rational_at(_q(U), S, q0)
    target, _ = uqsl2(S, target_q; max_terms=U.max_terms)
    return SpecializationMap(U, target, deepcopy(q0))
end

function (map::SpecializationMap)(x::UqSl2Elem)
    parent(x) === map.source || throw(ParentMismatchError())
    U = map.target
    C = elem_type(base_ring(U))
    terms = Dict{PBWMonomial,C}()
    for (m, coeff) in _terms(x)
        value = _evaluate_rational_at(coeff, base_ring(U), map.point)
        _add!(terms, m, value)
    end
    return UqSl2Elem(U, terms)
end
