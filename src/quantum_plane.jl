# Convention: X*Y = q*Y*X; the basis is Y^a*X^b.
struct PlaneMonomial
    a::Int
    b::Int
    function PlaneMonomial(a::Int, b::Int)
        a < 0 && throw(DomainError(a, "Y exponent must be nonnegative"))
        b < 0 && throw(DomainError(b, "X exponent must be nonnegative"))
        new(a, b)
    end
end

mutable struct QuantumPlaneAlgebra{C,R<:AA.Field} <: AA.NCRing
    const base::R
    const q::C
    const max_terms::Union{Nothing,Int}
end
Base.getproperty(P::QuantumPlaneAlgebra, name::Symbol) =
    name === :q ? deepcopy(getfield(P, :q)) : getfield(P, name)
_plane_q(P::QuantumPlaneAlgebra) = getfield(P, :q)

struct QuantumPlaneElem{C,R<:AA.Field} <: AA.NCRingElem
    parent::QuantumPlaneAlgebra{C,R}
    terms::Dict{PlaneMonomial,C}
    function QuantumPlaneElem(P::QuantumPlaneAlgebra{C,R}, terms::Dict{PlaneMonomial,C}) where {C,R}
        snapshot = Dict{PlaneMonomial,C}()
        for (m, value) in terms
            AA.parent(value) === base_ring(P) || throw(ArgumentError("coefficient belongs to another field"))
            !iszero(value) && (snapshot[m] = deepcopy(value))
        end
        _check_limit(P, length(snapshot))
        return new{C,R}(P, snapshot)
    end
end
_plane_terms(x::QuantumPlaneElem) = getfield(x, :terms)
Base.getproperty(x::QuantumPlaneElem, name::Symbol) =
    name === :terms ? throw(ArgumentError("use pbw_terms to read coefficients")) : getfield(x, name)
Base.propertynames(::QuantumPlaneElem, private::Bool=false) = private ? (:parent, :terms) : (:parent,)

"""
    quantum_plane(R, parameter; max_terms=nothing)

Construct the quantum plane `X*Y = q*Y*X` over the exact characteristic-zero
field `R`, with `q = R(parameter)`. Return the parent and `(X, Y)`; the PBW
basis is `Y^a X^b`. `q = 1` is treated as the commutative case and `q = 0` is
rejected. `max_terms` bounds the terms kept by each operation.
"""
function quantum_plane(R::AA.Field, parameter; max_terms::Union{Nothing,Integer}=nothing)
    is_exact_type(elem_type(R)) || throw(ArgumentError("coefficient field must be exact"))
    iszero(characteristic(R)) || throw(ArgumentError("coefficient field must have characteristic zero"))
    limit = if max_terms === nothing
        nothing
    else
        max_terms > 0 || throw(ArgumentError("max_terms must be positive"))
        _as_int(max_terms)
    end
    q = try
        R(parameter)
    catch
        throw(ArgumentError("parameter cannot be coerced into the coefficient field"))
    end
    iszero(q) && throw(DomainError(q, "q must be nonzero"))
    P = QuantumPlaneAlgebra{elem_type(R),typeof(R)}(R, deepcopy(q), limit)
    return P, gens(P)
end
quantum_plane(R, parameter; kwargs...) =
    throw(ArgumentError("coefficient parent must be an AbstractAlgebra field"))

base_ring(P::QuantumPlaneAlgebra) = P.base
base_ring(x::QuantumPlaneElem) = base_ring(parent(x))
parent(x::QuantumPlaneElem) = x.parent
elem_type(::Type{QuantumPlaneAlgebra{C,R}}) where {C,R} = QuantumPlaneElem{C,R}
elem_type(P::QuantumPlaneAlgebra) = elem_type(typeof(P))
parent_type(::Type{QuantumPlaneElem{C,R}}) where {C,R} = QuantumPlaneAlgebra{C,R}
base_ring_type(::Type{QuantumPlaneAlgebra{C,R}}) where {C,R} = R
characteristic(P::QuantumPlaneAlgebra) = characteristic(base_ring(P))
is_exact_type(::Type{QuantumPlaneElem{C,R}}) where {C,R} = is_exact_type(C)
quantum_parameter(P::QuantumPlaneAlgebra) = deepcopy(_plane_q(P))

function _check_limit(P::QuantumPlaneAlgebra, count::Int)
    limit = getfield(P, :max_terms)
    limit !== nothing && count > limit && throw(ComputationLimitError())
    return nothing
end

_plane_empty(P::QuantumPlaneAlgebra{C,R}) where {C,R} = QuantumPlaneElem(P, Dict{PlaneMonomial,C}())
zero(P::QuantumPlaneAlgebra) = _plane_empty(P)
one(P::QuantumPlaneAlgebra) = P(1)
zero(x::QuantumPlaneElem) = zero(parent(x))
one(x::QuantumPlaneElem) = one(parent(x))
iszero(x::QuantumPlaneElem) = isempty(_plane_terms(x))
isone(x::QuantumPlaneElem) = length(_plane_terms(x)) == 1 &&
    get(_plane_terms(x), PlaneMonomial(0,0), nothing) == one(base_ring(x))

function (P::QuantumPlaneAlgebra{C,R})(s) where {C,R}
    if s isa QuantumPlaneElem
        parent(s) === P || throw(ParentMismatchError())
        return copy(s)
    end
    value = try
        base_ring(P)(s)
    catch
        throw(ArgumentError("coefficient cannot be coerced into the base field"))
    end
    value isa C || throw(ArgumentError("coefficient has the wrong element type"))
    iszero(value) && return zero(P)
    return QuantumPlaneElem(P, Dict(PlaneMonomial(0,0) => value))
end

"""
    basis_monomial(P, a, b)

Return the basis element `Y^a X^b` of the quantum plane `P` (`a, b >= 0`).
"""
function basis_monomial(P::QuantumPlaneAlgebra, a::Integer, b::Integer)
    m = PlaneMonomial(_as_int(a), _as_int(b))
    return QuantumPlaneElem(P, Dict(m => one(base_ring(P))))
end
gens(P::QuantumPlaneAlgebra) = (basis_monomial(P,0,1), basis_monomial(P,1,0))
ngens(::QuantumPlaneAlgebra) = 2
"""
    pbw_coefficient(x, a, b)

Return the coefficient of `Y^a X^b` in `x`; a missing term gives the zero of
the base ring.
"""
pbw_coefficient(x::QuantumPlaneElem, a::Integer, b::Integer) =
    deepcopy(get(_plane_terms(x), PlaneMonomial(_as_int(a),_as_int(b)), zero(base_ring(x))))
_plane_order(m::PlaneMonomial) = (Base.Checked.checked_add(m.a,m.b), m.a, m.b)
"""
    pbw_terms(x)

Return the nonzero terms of `x` as a vector of `(PlaneMonomial, coefficient)`
pairs, sorted by `(a+b, a, b)` in descending order. Coefficients are copies.
"""
function pbw_terms(x::QuantumPlaneElem)
    keys_sorted = _sorted_keys(_plane_terms(x), _plane_order)
    return [(m, deepcopy(_plane_terms(x)[m])) for m in keys_sorted]
end
"""
    normal_form(x)

Return a copy of `x`. Elements are kept in PBW normal form at all times, so
this never changes the value.
"""
normal_form(x::QuantumPlaneElem) = copy(x)
copy(x::QuantumPlaneElem) = QuantumPlaneElem(parent(x),
    Dict(m => deepcopy(v) for (m,v) in _plane_terms(x)))
deepcopy_internal(x::QuantumPlaneElem, ::IdDict) = copy(x)

function _same_plane_parent(x::QuantumPlaneElem, y::QuantumPlaneElem)
    parent(x) === parent(y) || throw(ParentMismatchError())
    return parent(x)
end
function +(x::QuantumPlaneElem, y::QuantumPlaneElem)
    P = _same_plane_parent(x,y)
    terms = Dict(m => deepcopy(v) for (m,v) in _plane_terms(x))
    for (m,v) in _plane_terms(y)
        _add!(terms,m,v)
    end
    return QuantumPlaneElem(P, terms)
end
+(x::QuantumPlaneElem) = copy(x)
-(x::QuantumPlaneElem) = QuantumPlaneElem(parent(x),
    Dict(m => -v for (m,v) in _plane_terms(x)))
-(x::QuantumPlaneElem, y::QuantumPlaneElem) = x + (-y)

function _plane_scale(x::QuantumPlaneElem, s)
    P = parent(x)
    coeff = base_ring(P)(s)
    iszero(coeff) && return zero(P)
    return QuantumPlaneElem(P, Dict(m => coeff*v for (m,v) in _plane_terms(x)))
end
*(x::QuantumPlaneElem, s::Coefficient) = _plane_scale(x,s)
*(s::Coefficient, x::QuantumPlaneElem) = _plane_scale(x,s)
+(x::QuantumPlaneElem, s::Coefficient) = x + parent(x)(s)
+(s::Coefficient, x::QuantumPlaneElem) = x + s
-(x::QuantumPlaneElem, s::Coefficient) = x - parent(x)(s)
-(s::Coefficient, x::QuantumPlaneElem) = parent(x)(s) - x
function _plane_divide_scalar(x::QuantumPlaneElem, s)
    coeff = base_ring(x)(s)
    iszero(coeff) && throw(DivideError())
    return _plane_scale(x, inv(coeff))
end
/(x::QuantumPlaneElem, s::Integer) = _plane_divide_scalar(x,s)
/(x::QuantumPlaneElem, s::Rational) = _plane_divide_scalar(x,s)
/(x::QuantumPlaneElem, s::AA.RingElem) = _plane_divide_scalar(x,s)

function *(x::QuantumPlaneElem, y::QuantumPlaneElem)
    P = _same_plane_parent(x,y)
    terms = Dict{PlaneMonomial,elem_type(base_ring(P))}()
    for (left,lv) in _plane_terms(x), (right,rv) in _plane_terms(y)
        m = PlaneMonomial(Base.Checked.checked_add(left.a,right.a),
                          Base.Checked.checked_add(left.b,right.b))
        exponent = Base.Checked.checked_mul(left.b,right.a)
        _add!(terms, m, lv*rv*_plane_q(P)^exponent)
    end
    return QuantumPlaneElem(P, terms)
end
function inv(x::QuantumPlaneElem)
    length(_plane_terms(x)) == 1 || throw(UnsupportedInverseError())
    m, value = first(_plane_terms(x))
    m == PlaneMonomial(0,0) || throw(UnsupportedInverseError())
    return parent(x)(inv(value))
end
function ^(x::QuantumPlaneElem, n::Integer)
    if n < 0
        magnitude = Base.Checked.checked_neg(n)
        return inv(x)^magnitude
    end
    result = one(x)
    factor = x
    while n > 0
        isodd(n) && (result = result*factor)
        n = n >> 1
        n > 0 && (factor = factor*factor)
    end
    return result
end

==(x::QuantumPlaneElem, y::QuantumPlaneElem) = parent(x) === parent(y) && _plane_terms(x) == _plane_terms(y)
isequal(x::QuantumPlaneElem, y::QuantumPlaneElem) = parent(x) === parent(y) && isequal(_plane_terms(x),_plane_terms(y))
isequal(::QuantumPlaneElem, ::Coefficient) = false
isequal(::Coefficient, ::QuantumPlaneElem) = false
==(x::QuantumPlaneElem, s::Coefficient) = x == parent(x)(s)
==(s::Coefficient, x::QuantumPlaneElem) = x == s
function hash(x::QuantumPlaneElem, h::UInt)
    h = hash(objectid(parent(x)),hash(:QuantumPlaneElem,h))
    for (m,value) in pbw_terms(x)
        h = hash((m.a,m.b,value),h)
    end
    return h
end
"""
    commutator(x, y)

Return `x*y - y*x`.
"""
commutator(x::QuantumPlaneElem,y::QuantumPlaneElem) = x*y-y*x

function show(io::IO, P::QuantumPlaneAlgebra)
    print(io,"Quantum plane over ")
    show(io,base_ring(P))
    print(io," with q = ")
    show(io,_plane_q(P))
end
function show(io::IO, x::QuantumPlaneElem)
    iszero(x) && return print(io,"0")
    for (i,(m,value)) in enumerate(pbw_terms(x))
        i > 1 && print(io," + ")
        factors = String[]
        m.a > 0 && push!(factors,m.a == 1 ? "Y" : "Y^$(m.a)")
        m.b > 0 && push!(factors,m.b == 1 ? "X" : "X^$(m.b)")
        if isempty(factors) || !isone(value)
            show(io,value)
            !isempty(factors) && print(io,"*")
        end
        print(io,join(factors,"*"))
    end
end
