# A free-word normalizer used only in tests. It never calls the package's
# multiplication or its PBW right-multiplication recurrence.
const Word = Tuple{Vararg{Symbol}}
const REDEXES = Dict(
    (:K, :F) => ((:F, :K),),
    (:L, :F) => ((:F, :L),),
    (:E, :K) => ((:K, :E),),
    (:E, :L) => ((:L, :E),),
    (:E, :F) => ((:F, :E), (:K,), (:L,)),
    (:K, :L) => ((),),
    (:L, :K) => ((),),
)

function reference_normalize(word::Word, q; choose_right=false)
    C = typeof(q)
    terms = Dict{Word,C}(word => one(q))
    while true
        changed = false
        next_terms = Dict{Word,C}()
        for (w, coeff) in terms
            positions = [i for i in 1:(length(w)-1) if haskey(REDEXES, (w[i], w[i+1]))]
            if isempty(positions)
                next_terms[w] = get(next_terms, w, zero(q)) + coeff
                continue
            end
            changed = true
            i = choose_right ? last(positions) : first(positions)
            pair = (w[i], w[i+1])
            replacements = REDEXES[pair]
            multipliers = pair == (:K, :F) ? (q^-2,) :
                pair == (:L, :F) ? (q^2,) :
                pair == (:E, :K) ? (q^-2,) :
                pair == (:E, :L) ? (q^2,) :
                pair == (:E, :F) ? (one(q), inv(q-inv(q)), -inv(q-inv(q))) :
                (one(q),)
            for (replacement, multiplier) in zip(replacements, multipliers)
                key = (w[1:i-1]..., replacement..., w[i+2:end]...)
                next_terms[key] = get(next_terms, key, zero(q)) + coeff*multiplier
            end
        end
        filter!(pair -> !iszero(last(pair)), next_terms)
        terms = next_terms
        changed || return terms
    end
end

function reference_as_pbw(U, terms)
    result = zero(U)
    for (word, coeff) in terms
        a = count(==(:F), word)
        b = count(==(:K), word) - count(==(:L), word)
        c = count(==(:E), word)
        result += coeff*basis_monomial(U, a, b, c)
    end
    return result
end

function check_reference_rewriter(U, q)
    E, F, K, Ki = AA.gens(U)
    generators = Dict(:E => E, :F => F, :K => K, :L => Ki)
    overlaps = [(:E,:K,:F), (:E,:L,:F), (:E,:K,:L), (:E,:L,:K),
        (:K,:L,:F), (:L,:K,:F), (:K,:L,:K), (:L,:K,:L)]
    for word in overlaps
        @test reference_normalize(word, q) == reference_normalize(word, q; choose_right=true)
    end
    for len in 0:4, chars in Iterators.product(ntuple(_ -> (:E,:F,:K,:L), len)...)
        word = Tuple(chars)
        reference = reference_as_pbw(U, reference_normalize(word, q))
        actual = foldl(*, (generators[c] for c in word); init=one(U))
        @test reference == actual
    end
end
