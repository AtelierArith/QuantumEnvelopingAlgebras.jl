import AbstractAlgebra as AA
using QuantumEnvelopingAlgebras

R, q = AA.rational_function_field(AA.QQ, :q)
U, (E, F, K, Ki) = uqsl2(R, q)

cases = [
    ("E^2 F^2", () -> E^2 * F^2),
    ("E^4 F^4", () -> E^4 * F^4),
    ("(E+F)^6", () -> (E+F)^6),
    ("negative K power", () -> E*K^(-12)*F),
    ("sparse product", () -> (E^2+K^(-3)+F^2)*(F^2+K^4+E^2)),
]

println("Julia ", VERSION)
println("case,time_s,bytes,terms")
for (name, run) in cases
    run() # compile before measuring
    result = @timed run()
    println(name, ",", result.time, ",", result.bytes, ",", length(pbw_terms(result.value)))
end

P, (X, Y) = quantum_plane(R, q)
plane_cases = [
    ("plane (X+Y)^8", () -> (X+Y)^8),
    ("plane sparse product", () -> (X^4+Y^4+one(P))*(Y^4+X^4+one(P))),
]
for (name, run) in plane_cases
    run()
    result = @timed run()
    println(name, ",", result.time, ",", result.bytes, ",", length(pbw_terms(result.value)))
end
