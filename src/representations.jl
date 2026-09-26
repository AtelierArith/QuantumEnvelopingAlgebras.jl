"A type-1 highest-weight representation with basis v_0,...,v_n."
struct HighestWeightRepresentation{U,M,A}
    source::U
    space::M
    E::A
    F::A
    K::A
    Ki::A
end

Base.getproperty(rho::HighestWeightRepresentation, name::Symbol) =
    name in (:E, :F, :K, :Ki) ? deepcopy(getfield(rho, name)) : getfield(rho, name)
_matrix(rho::HighestWeightRepresentation, name::Symbol) = getfield(rho, name)

"""
    highest_weight_representation(U, n)

Return the type-1 highest-weight representation of dimension `n+1` as an
`AbstractAlgebra` matrix representation. Call `rho(x)` to map an element of `U`
to a matrix; the basis matrices are available as `rho.E`, `rho.F`, `rho.K` and
`rho.Ki`.
"""
function highest_weight_representation(U::UqSl2Algebra, n::Integer)
    n >= 0 || throw(DomainError(n, "highest weight must be nonnegative"))
    weight = _as_int(n)
    dimension = Base.Checked.checked_add(weight, 1)
    M = AA.matrix_ring(base_ring(U), dimension)
    E, F, K, Ki = M(), M(), M(), M()
    q = _q(U)
    for j in 0:weight
        K[j+1,j+1] = q^(weight - 2j)
        Ki[j+1,j+1] = q^(-weight + 2j)
        j < weight && (F[j+2,j+1] = one(base_ring(U)))
        j > 0 && (E[j,j+1] = qinteger(q,j)*qinteger(q,weight-j+1))
    end
    return HighestWeightRepresentation(U, M, E, F, K, Ki)
end

function (rho::HighestWeightRepresentation)(x::UqSl2Elem)
    parent(x) === rho.source || throw(ParentMismatchError())
    result = zero(rho.space)
    for (m, coeff) in _terms(x)
        Kpower = m.b >= 0 ? _matrix(rho, :K)^m.b : _matrix(rho, :Ki)^Base.Checked.checked_neg(m.b)
        result += coeff*(_matrix(rho, :F)^m.a * Kpower * _matrix(rho, :E)^m.c)
    end
    return result
end
