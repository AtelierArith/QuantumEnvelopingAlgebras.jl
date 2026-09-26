module QuantumEnvelopingAlgebras

import AbstractAlgebra as AA
import AbstractAlgebra: parent, base_ring, elem_type, parent_type, base_ring_type,
    characteristic, is_exact_type, gens, ngens
import Base: zero, one, iszero, isone, +, -, *, /, ^, inv, ==, isequal, hash,
    copy, deepcopy_internal, show

export uqsl2, PBWMonomial, quantum_parameter,
    basis_monomial, pbw_coefficient, pbw_terms, normal_form, commutator,
    qinteger, casimir, ParentMismatchError, UnsupportedInverseError,
    ComputationLimitError, tensor_power, tensor, coproduct_map,
    coproduct_on_slot, counit, antipode, highest_weight_representation,
    specialization_map, SpecializationError, quantum_plane

"Raised when an operation mixes elements of different parents."
struct ParentMismatchError <: Exception end
show(io::IO, ::ParentMismatchError) = print(io, "elements have different U_q(sl_2) parents")
"Raised when an inverse is requested for an element that is not `s*K^b`."
struct UnsupportedInverseError <: Exception end
show(io::IO, ::UnsupportedInverseError) = print(io, "inverse is not supported for this element")
"Raised when an operation would exceed the parent's `max_terms` budget."
struct ComputationLimitError <: Exception end
show(io::IO, ::ComputationLimitError) = print(io, "computation limit exceeded")
"Raised when a specialization's coefficient denominator vanishes at the point."
struct SpecializationError <: Exception
    message::String
end
Base.showerror(io::IO, error::SpecializationError) = print(io, error.message)

"""
    PBWMonomial(a, b, c)

PBW basis key `F^a K^b E^c` with `a, c >= 0` and arbitrary `b`. This is the
monomial type in the pairs returned by [`pbw_terms`](@ref).
"""
struct PBWMonomial
    a::Int
    b::Int
    c::Int
    function PBWMonomial(a::Int, b::Int, c::Int)
        a < 0 && throw(DomainError(a, "F exponent must be nonnegative"))
        c < 0 && throw(DomainError(c, "E exponent must be nonnegative"))
        new(a, b, c)
    end
end

mutable struct UqSl2Algebra{C,R<:AA.Field} <: AA.NCRing
    const base::R
    const q::C
    const qinv::C
    const dinv::C
    const max_terms::Union{Nothing,Int}
end

Base.getproperty(U::UqSl2Algebra, name::Symbol) =
    name in (:q, :qinv, :dinv) ? deepcopy(getfield(U, name)) : getfield(U, name)
_q(U::UqSl2Algebra) = getfield(U, :q)
_qinv(U::UqSl2Algebra) = getfield(U, :qinv)
_dinv(U::UqSl2Algebra) = getfield(U, :dinv)

struct UqSl2Elem{C,R<:AA.Field} <: AA.NCRingElem
    parent::UqSl2Algebra{C,R}
    terms::Dict{PBWMonomial,C}
    function UqSl2Elem(U::UqSl2Algebra{C,R}, terms::Dict{PBWMonomial,C}) where {C,R}
        snapshot = Dict{PBWMonomial,C}()
        for (m, value) in terms
            AA.parent(value) === base_ring(U) || throw(ArgumentError("coefficient belongs to another field"))
            !iszero(value) && (snapshot[m] = deepcopy(value))
        end
        _check_limit(U, length(snapshot))
        return new{C,R}(U, snapshot)
    end
end

_terms(x::UqSl2Elem) = getfield(x, :terms)
Base.getproperty(x::UqSl2Elem, name::Symbol) =
    name === :terms ? throw(ArgumentError("use pbw_terms to read coefficients")) : getfield(x, name)
Base.propertynames(::UqSl2Elem, private::Bool=false) = private ? (:parent, :terms) : (:parent,)

"""
    uqsl2(R::AbstractAlgebra.Field, parameter; max_terms=nothing)

Construct the Drinfeld–Jimbo algebra `U_q(sl_2)` over the exact
characteristic-zero field `R`, with `q = R(parameter)`.

Return the parent and the generators `(E, F, K, K_i)`. Each call creates a new
parent, and elements of different parents do not mix. `max_terms` caps the
number of PBW terms retained by each operation; exceeding it throws
[`ComputationLimitError`](@ref).
"""
function uqsl2(R::AA.Field, parameter; max_terms::Union{Nothing,Integer}=nothing)
    limit = if max_terms === nothing
        nothing
    else
        max_terms > 0 || throw(ArgumentError("max_terms must be positive"))
        _as_int(max_terms)
    end
    is_exact_type(elem_type(R)) || throw(ArgumentError("coefficient field must be exact"))
    iszero(characteristic(R)) || throw(ArgumentError("coefficient field must have characteristic zero"))
    q = try
        R(parameter)
    catch
        throw(ArgumentError("parameter cannot be coerced into the coefficient field"))
    end
    iszero(q) && throw(DomainError(q, "q must be nonzero"))
    qinv = inv(q)
    d = q - qinv
    iszero(d) && throw(DomainError(q, "q - q^(-1) must be nonzero"))
    U = UqSl2Algebra{elem_type(R),typeof(R)}(R, deepcopy(q), deepcopy(qinv), deepcopy(inv(d)), limit)
    return U, gens(U)
end
uqsl2(R, parameter; kwargs...) = throw(ArgumentError("coefficient parent must be an AbstractAlgebra field"))

base_ring(U::UqSl2Algebra) = U.base
base_ring(x::UqSl2Elem) = base_ring(parent(x))
parent(x::UqSl2Elem) = x.parent
elem_type(::Type{UqSl2Algebra{C,R}}) where {C,R} = UqSl2Elem{C,R}
elem_type(U::UqSl2Algebra) = elem_type(typeof(U))
parent_type(::Type{UqSl2Elem{C,R}}) where {C,R} = UqSl2Algebra{C,R}
base_ring_type(::Type{UqSl2Algebra{C,R}}) where {C,R} = R
characteristic(U::UqSl2Algebra) = characteristic(base_ring(U))
is_exact_type(::Type{UqSl2Elem{C,R}}) where {C,R} = is_exact_type(C)
"""
    quantum_parameter(U)

Return a copy of the deformation parameter `q` of the parent `U`.
"""
quantum_parameter(U::UqSl2Algebra) = deepcopy(_q(U))
function _check_limit(U::UqSl2Algebra, count::Int)
    limit = getfield(U, :max_terms)
    limit !== nothing && count > limit && throw(ComputationLimitError())
    return nothing
end
function _checked(x::UqSl2Elem)
    _check_limit(parent(x), length(_terms(x)))
    return x
end

_empty(U::UqSl2Algebra{C,R}) where {C,R} = UqSl2Elem(U, Dict{PBWMonomial,C}())
function _add!(terms::Dict{M,C}, m::M, value) where {M,C}
    iszero(value) && return terms
    newvalue = haskey(terms, m) ? terms[m] + value : value
    if iszero(newvalue)
        delete!(terms, m)
    else
        terms[m] = deepcopy(newvalue)
    end
    return terms
end
function _sorted_keys(terms, order)
    keys_sorted = sort!(collect(keys(terms)); by=order, rev=true)
    foreach(order, keys_sorted) # also checks overflow for a single key
    return keys_sorted
end

function (U::UqSl2Algebra{C,R})(s) where {C,R}
    if s isa UqSl2Elem
        parent(s) === U || throw(ParentMismatchError())
        return copy(s)
    end
    value = try
        U.base(s)
    catch
        throw(ArgumentError("coefficient cannot be coerced into the base field"))
    end
    value isa C || throw(ArgumentError("coefficient has the wrong element type"))
    iszero(value) && return _empty(U)
    return UqSl2Elem(U, Dict(PBWMonomial(0, 0, 0) => deepcopy(value)))
end

zero(U::UqSl2Algebra) = _empty(U)
one(U::UqSl2Algebra) = U(1)
zero(x::UqSl2Elem) = zero(parent(x))
one(x::UqSl2Elem) = one(parent(x))
iszero(x::UqSl2Elem) = isempty(_terms(x))
isone(x::UqSl2Elem) = length(_terms(x)) == 1 &&
    get(_terms(x), PBWMonomial(0, 0, 0), nothing) == one(base_ring(x))

function _as_int(n::Integer)
    (n < typemin(Int) || n > typemax(Int)) && throw(OverflowError("PBW exponent exceeds Int range"))
    return Int(n)
end
"""
    basis_monomial(U, a, b, c)

Return the basis element `F^a K^b E^c` of `U` (`a, c >= 0`).
"""
function basis_monomial(U::UqSl2Algebra, a::Integer, b::Integer, c::Integer)
    m = PBWMonomial(_as_int(a), _as_int(b), _as_int(c))
    return UqSl2Elem(U, Dict(m => one(base_ring(U))))
end
gens(U::UqSl2Algebra) = (basis_monomial(U, 0, 0, 1), basis_monomial(U, 1, 0, 0),
    basis_monomial(U, 0, 1, 0), basis_monomial(U, 0, -1, 0))
ngens(::UqSl2Algebra) = 4
"""
    pbw_coefficient(x, a, b, c)

Return the coefficient of `F^a K^b E^c` in `x`; a missing term gives the zero
of the base ring.
"""
pbw_coefficient(x::UqSl2Elem, a::Integer, b::Integer, c::Integer) =
    deepcopy(get(_terms(x), PBWMonomial(_as_int(a), _as_int(b), _as_int(c)), zero(base_ring(x))))
_order(m::PBWMonomial) = (Base.Checked.checked_add(m.a, m.c), m.a, m.b, m.c)
"""
    pbw_terms(x)

Return the nonzero terms of `x` as a vector of `(PBWMonomial, coefficient)`
pairs, sorted by `(a+c, a, b, c)` in descending order. Coefficients are copies.
"""
function pbw_terms(x::UqSl2Elem)
    keys_sorted = _sorted_keys(_terms(x), _order)
    return [(m, deepcopy(_terms(x)[m])) for m in keys_sorted]
end
"""
    normal_form(x)

Return a copy of `x`. Elements are kept in PBW normal form at all times, so
this never changes the value.
"""
normal_form(x::UqSl2Elem) = copy(x)
copy(x::UqSl2Elem) = UqSl2Elem(parent(x), Dict(m => deepcopy(v) for (m, v) in _terms(x)))
deepcopy_internal(x::UqSl2Elem, ::IdDict) = copy(x)

function _same_parent(x::UqSl2Elem, y::UqSl2Elem)
    parent(x) === parent(y) || throw(ParentMismatchError())
    return parent(x)
end
function +(x::UqSl2Elem, y::UqSl2Elem)
    _same_parent(x, y)
    result = copy(x)
    for (m, v) in _terms(y)
        _add!(_terms(result), m, v)
    end
    return _checked(result)
end
+(x::UqSl2Elem) = copy(x)
function -(x::UqSl2Elem)
    result = zero(x)
    for (m, v) in _terms(x)
        _add!(_terms(result), m, -v)
    end
    return _checked(result)
end
-(x::UqSl2Elem, y::UqSl2Elem) = x + (-y)

function _scale(x::UqSl2Elem, s)
    value = base_ring(x)(s)
    result = zero(x)
    for (m, v) in _terms(x)
        _add!(_terms(result), m, value * v)
    end
    return _checked(result)
end
const Coefficient = Union{Integer,Rational,AA.RingElem}
+(x::UqSl2Elem, s::Coefficient) = x + parent(x)(s)
+(s::Coefficient, x::UqSl2Elem) = x + s
-(x::UqSl2Elem, s::Coefficient) = x - parent(x)(s)
-(s::Coefficient, x::UqSl2Elem) = parent(x)(s) - x
*(x::UqSl2Elem, s::Coefficient) = _scale(x, s)
*(s::Coefficient, x::UqSl2Elem) = _scale(x, s)
function _divide_scalar(x::UqSl2Elem, s)
    value = base_ring(x)(s)
    iszero(value) && throw(DivideError())
    return _scale(x, inv(value))
end
/(x::UqSl2Elem, s::Integer) = _divide_scalar(x, s)
/(x::UqSl2Elem, s::Rational) = _divide_scalar(x, s)
/(x::UqSl2Elem, s::AA.RingElem) = _divide_scalar(x, s)

function _right_F(x::UqSl2Elem)
    U = parent(x)
    result = zero(U)
    for (m, v) in _terms(x)
        _add!(_terms(result), PBWMonomial(Base.Checked.checked_add(m.a, 1), m.b, m.c),
            v * _q(U)^Base.Checked.checked_mul(-2, m.b))
        if m.c > 0
            factor = v * qinteger(_q(U), m.c) * _dinv(U)
            _add!(_terms(result), PBWMonomial(m.a, Base.Checked.checked_add(m.b, 1), m.c - 1),
                factor * _q(U)^(1 - m.c))
            _add!(_terms(result), PBWMonomial(m.a, Base.Checked.checked_sub(m.b, 1), m.c - 1),
                -factor * _q(U)^(m.c - 1))
        end
    end
    return _checked(result)
end
function _right_Kpow(x::UqSl2Elem, b::Int)
    U = parent(x)
    result = zero(U)
    for (m, v) in _terms(x)
        exponent = Base.Checked.checked_mul(-2, Base.Checked.checked_mul(m.c, b))
        _add!(_terms(result), PBWMonomial(m.a, Base.Checked.checked_add(m.b, b), m.c),
            v * _q(U)^exponent)
    end
    return _checked(result)
end
function _right_Epow(x::UqSl2Elem, c::Int)
    result = zero(x)
    for (m, v) in _terms(x)
        _add!(_terms(result), PBWMonomial(m.a, m.b, Base.Checked.checked_add(m.c, c)), v)
    end
    return _checked(result)
end
function *(x::UqSl2Elem{C}, y::UqSl2Elem) where C
    U = _same_parent(x, y)
    terms = Dict{PBWMonomial,C}()
    for (left, lv) in _terms(x), (right, rv) in _terms(y)
        piece = basis_monomial(U, left.a, left.b, left.c)
        for _ in 1:right.a
            piece = _right_F(piece)
        end
        piece = _right_Epow(_right_Kpow(piece, right.b), right.c)
        scalar = lv*rv
        for (m, value) in _terms(piece)
            _add!(terms, m, scalar*value)
        end
        _check_limit(U, length(terms))
    end
    return UqSl2Elem(U, terms)
end
"""
    qinteger(q, n)

Return the quantum integer `[n]_q = sum_{j=0}^{n-1} q^(n-1-2j)` for `n >= 0`.
"""
function qinteger(q, n::Integer)
    n < 0 && throw(DomainError(n, "q integer index must be nonnegative"))
    result = zero(q)
    for j in 0:(n - 1)
        result += q^(n - 1 - 2j)
    end
    return result
end
function inv(x::UqSl2Elem)
    length(_terms(x)) == 1 || throw(UnsupportedInverseError())
    m, s = first(_terms(x))
    (m.a == 0 && m.c == 0 && !iszero(s)) || throw(UnsupportedInverseError())
    return _scale(basis_monomial(parent(x), 0, Base.Checked.checked_neg(m.b), 0), inv(s))
end
function ^(x::UqSl2Elem, n::Integer)
    if n < 0
        magnitude = Base.Checked.checked_neg(n)
        return inv(x)^magnitude
    end
    result = one(x)
    factor = x
    while n > 0
        isodd(n) && (result = result * factor)
        n = n >> 1
        n > 0 && (factor = factor * factor)
    end
    return result
end
==(x::UqSl2Elem, y::UqSl2Elem) = parent(x) === parent(y) && _terms(x) == _terms(y)
isequal(x::UqSl2Elem, y::UqSl2Elem) = parent(x) === parent(y) && isequal(_terms(x), _terms(y))
isequal(::UqSl2Elem, ::Coefficient) = false
isequal(::Coefficient, ::UqSl2Elem) = false
==(x::UqSl2Elem, s::Coefficient) = x == parent(x)(s)
==(s::Coefficient, x::UqSl2Elem) = x == s
function hash(x::UqSl2Elem, h::UInt)
    h = hash(objectid(parent(x)), hash(:UqSl2Elem, h))
    for (m, v) in pbw_terms(x)
        h = hash((m.a, m.b, m.c, v), h)
    end
    return h
end
"""
    commutator(x, y)

Return `x*y - y*x`.
"""
commutator(x::UqSl2Elem, y::UqSl2Elem) = x*y - y*x
"""
    casimir(U)

Return the Casimir element `F*E + (q*K + q^-1*K^-1)/(q - q^-1)^2`, which is
central in `U`.
"""
function casimir(U::UqSl2Algebra)
    E, F, K, Ki = gens(U)
    return F*E + (_q(U)*K + _qinv(U)*Ki) * _dinv(U)^2
end
function show(io::IO, U::UqSl2Algebra)
    print(io, "U_q(sl_2) over ")
    show(io, base_ring(U))
    print(io, " with q = ")
    show(io, _q(U))
end
function show(io::IO, x::UqSl2Elem)
    iszero(x) && return print(io, "0")
    for (i, (m, v)) in enumerate(pbw_terms(x))
        i > 1 && print(io, " + ")
        factors = String[]
        m.a > 0 && push!(factors, m.a == 1 ? "F" : "F^$(m.a)")
        m.b != 0 && push!(factors, m.b == 1 ? "K" : "K^$(m.b)")
        m.c > 0 && push!(factors, m.c == 1 ? "E" : "E^$(m.c)")
        if isempty(factors) || !isone(v)
            show(io, v)
            !isempty(factors) && print(io, "*")
        end
        print(io, join(factors, "*"))
    end
end

include("tensor.jl")
include("hopf.jl")
include("representations.jl")
include("specialization.jl")
include("quantum_plane.jl")

end # module
