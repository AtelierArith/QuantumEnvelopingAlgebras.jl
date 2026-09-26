"A parent for the N-fold tensor power of one particular U_q(sl_2) parent."
mutable struct UqSl2TensorPower{N,C,R<:AA.Field} <: AA.NCRing
    const algebra::UqSl2Algebra{C,R}
end

struct UqSl2TensorElem{N,C,R<:AA.Field} <: AA.NCRingElem
    parent::UqSl2TensorPower{N,C,R}
    terms::Dict{NTuple{N,PBWMonomial},C}
    function UqSl2TensorElem(T::UqSl2TensorPower{N,C,R},
                             terms::Dict{NTuple{N,PBWMonomial},C}) where {N,C,R}
        snapshot = Dict{NTuple{N,PBWMonomial},C}()
        for (key, value) in terms
            AA.parent(value) === base_ring(T) || throw(ArgumentError("coefficient belongs to another field"))
            !iszero(value) && (snapshot[key] = deepcopy(value))
        end
        _check_limit(T.algebra, length(snapshot))
        return new{N,C,R}(T, snapshot)
    end
end

_tensor_terms(x::UqSl2TensorElem) = getfield(x, :terms)
Base.getproperty(x::UqSl2TensorElem, name::Symbol) =
    name === :terms ? throw(ArgumentError("use tensor_terms to read coefficients")) : getfield(x, name)
Base.propertynames(::UqSl2TensorElem, private::Bool=false) = private ? (:parent, :terms) : (:parent,)

"""
    tensor_power(U, n)

Return the `n`-fold tensor power parent of the algebra `U` (`n > 0`). Each call
creates a new parent.
"""
function tensor_power(U::UqSl2Algebra{C,R}, n::Integer) where {C,R}
    n > 0 || throw(ArgumentError("tensor power must be positive"))
    N = _as_int(n)
    return UqSl2TensorPower{N,C,R}(U)
end

parent(x::UqSl2TensorElem) = x.parent
base_ring(T::UqSl2TensorPower) = base_ring(T.algebra)
base_ring(x::UqSl2TensorElem) = base_ring(parent(x))
elem_type(::Type{UqSl2TensorPower{N,C,R}}) where {N,C,R} = UqSl2TensorElem{N,C,R}
elem_type(T::UqSl2TensorPower) = elem_type(typeof(T))
parent_type(::Type{UqSl2TensorElem{N,C,R}}) where {N,C,R} = UqSl2TensorPower{N,C,R}
base_ring_type(::Type{UqSl2TensorPower{N,C,R}}) where {N,C,R} = R
characteristic(T::UqSl2TensorPower) = characteristic(T.algebra)
is_exact_type(::Type{UqSl2TensorElem{N,C,R}}) where {N,C,R} = is_exact_type(C)

_tensor_empty(T::UqSl2TensorPower{N,C,R}) where {N,C,R} =
    UqSl2TensorElem(T, Dict{NTuple{N,PBWMonomial},C}())
zero(T::UqSl2TensorPower) = _tensor_empty(T)
zero(x::UqSl2TensorElem) = zero(parent(x))
function one(T::UqSl2TensorPower{N,C,R}) where {N,C,R}
    key = ntuple(_ -> PBWMonomial(0, 0, 0), N)
    return UqSl2TensorElem(T, Dict(key => one(base_ring(T))))
end
one(x::UqSl2TensorElem) = one(parent(x))
iszero(x::UqSl2TensorElem) = isempty(_tensor_terms(x))
isone(x::UqSl2TensorElem) = x == one(parent(x))

function (T::UqSl2TensorPower)(s)
    if s isa UqSl2TensorElem
        parent(s) === T || throw(ParentMismatchError())
        return copy(s)
    end
    value = base_ring(T)(s)
    iszero(value) && return zero(T)
    key = ntuple(_ -> PBWMonomial(0, 0, 0), length(T))
    return UqSl2TensorElem(T, Dict(key => value))
end
Base.length(::UqSl2TensorPower{N}) where N = N

"""
    tensor(T, x1, ..., xn)

Embed elements `x1, ..., xn` of `T.algebra` as the pure tensor
`x1 ⊗ ... ⊗ xn` in the `n`-fold tensor power `T`.
"""
function tensor(T::UqSl2TensorPower{N,C,R}, xs::Vararg{UqSl2Elem,N}) where {N,C,R}
    all(parent(x) === T.algebra for x in xs) || throw(ParentMismatchError())
    terms = Dict{NTuple{N,PBWMonomial},C}()
    for parts in Iterators.product((collect(_terms(x)) for x in xs)...)
        key = ntuple(i -> first(parts[i]), N)
        coeff = foldl(*, (last(parts[i]) for i in 1:N); init=one(base_ring(T)))
        _add!(terms, key, coeff)
    end
    return UqSl2TensorElem(T, terms)
end

function tensor_terms(x::UqSl2TensorElem)
    keys_sorted = _sorted_keys(_tensor_terms(x), k -> map(_order, k))
    return [(k, deepcopy(_tensor_terms(x)[k])) for k in keys_sorted]
end
copy(x::UqSl2TensorElem) = UqSl2TensorElem(parent(x),
    Dict(k => deepcopy(v) for (k, v) in _tensor_terms(x)))
deepcopy_internal(x::UqSl2TensorElem, ::IdDict) = copy(x)

function _same_tensor_parent(x::UqSl2TensorElem, y::UqSl2TensorElem)
    parent(x) === parent(y) || throw(ParentMismatchError())
    return parent(x)
end
function +(x::UqSl2TensorElem, y::UqSl2TensorElem)
    T = _same_tensor_parent(x, y)
    terms = copy(_tensor_terms(x))
    for (key, value) in _tensor_terms(y)
        _add!(terms, key, value)
    end
    return UqSl2TensorElem(T, terms)
end
+(x::UqSl2TensorElem) = copy(x)
function -(x::UqSl2TensorElem)
    T = parent(x)
    return UqSl2TensorElem(T, Dict(k => -v for (k, v) in _tensor_terms(x)))
end
-(x::UqSl2TensorElem, y::UqSl2TensorElem) = x + (-y)

function _tensor_scale(x::UqSl2TensorElem, s)
    T = parent(x)
    coeff = base_ring(T)(s)
    iszero(coeff) && return zero(T)
    return UqSl2TensorElem(T, Dict(k => coeff*v for (k, v) in _tensor_terms(x)))
end
*(x::UqSl2TensorElem, s::Coefficient) = _tensor_scale(x, s)
*(s::Coefficient, x::UqSl2TensorElem) = _tensor_scale(x, s)
function _tensor_divide_scalar(x::UqSl2TensorElem, s)
    coeff = base_ring(x)(s)
    iszero(coeff) && throw(DivideError())
    return _tensor_scale(x, inv(coeff))
end
/(x::UqSl2TensorElem, s::Integer) = _tensor_divide_scalar(x, s)
/(x::UqSl2TensorElem, s::Rational) = _tensor_divide_scalar(x, s)
/(x::UqSl2TensorElem, s::AA.RingElem) = _tensor_divide_scalar(x, s)
+(x::UqSl2TensorElem, s::Coefficient) = x + parent(x)(s)
+(s::Coefficient, x::UqSl2TensorElem) = x + s
-(x::UqSl2TensorElem, s::Coefficient) = x - parent(x)(s)
-(s::Coefficient, x::UqSl2TensorElem) = parent(x)(s) - x

function *(x::UqSl2TensorElem, y::UqSl2TensorElem)
    T = _same_tensor_parent(x, y)
    U = T.algebra
    N = length(T)
    result = zero(T)
    for (left, lv) in _tensor_terms(x), (right, rv) in _tensor_terms(y)
        factors = ntuple(i -> basis_monomial(U, left[i].a, left[i].b, left[i].c) *
                              basis_monomial(U, right[i].a, right[i].b, right[i].c), N)
        result += (lv*rv)*tensor(T, factors...)
    end
    return result
end
function ^(x::UqSl2TensorElem, n::Integer)
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
function inv(x::UqSl2TensorElem{N}) where N
    length(_tensor_terms(x)) == 1 || throw(UnsupportedInverseError())
    key, coeff = first(_tensor_terms(x))
    all(m.a == 0 && m.c == 0 for m in key) || throw(UnsupportedInverseError())
    inverse_key = ntuple(i -> PBWMonomial(0, Base.Checked.checked_neg(key[i].b), 0), N)
    return UqSl2TensorElem(parent(x), Dict(inverse_key => inv(coeff)))
end

==(x::UqSl2TensorElem, y::UqSl2TensorElem) =
    parent(x) === parent(y) && _tensor_terms(x) == _tensor_terms(y)
isequal(x::UqSl2TensorElem, y::UqSl2TensorElem) =
    parent(x) === parent(y) && isequal(_tensor_terms(x), _tensor_terms(y))
function hash(x::UqSl2TensorElem, h::UInt)
    h = hash(objectid(parent(x)), hash(:UqSl2TensorElem, h))
    for (key, value) in tensor_terms(x)
        h = hash((map(m -> (m.a, m.b, m.c), key), value), h)
    end
    return h
end

function show(io::IO, T::UqSl2TensorPower{N}) where N
    print(io, "", N, "-fold tensor power of ")
    show(io, T.algebra)
end
function show(io::IO, x::UqSl2TensorElem)
    iszero(x) && return print(io, "0")
    for (i, (key, value)) in enumerate(tensor_terms(x))
        i > 1 && print(io, " + ")
        if !isone(value)
            show(io, value)
            print(io, "*")
        end
        print(io, "(")
        for (j, m) in enumerate(key)
            j > 1 && print(io, " ⊗ ")
            show(io, basis_monomial(parent(x).algebra, m.a, m.b, m.c))
        end
        print(io, ")")
    end
end
