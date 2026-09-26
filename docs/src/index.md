# QuantumEnvelopingAlgebras.jl

`U_q(sl_2)` と量子平面を厳密な PBW 正規形で計算する Julia パッケージです。係数体には
[AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl) を使います。

以下の例は `import AbstractAlgebra as AA` と `using QuantumEnvelopingAlgebras`
が済んでいるものとします。

```@meta
DocTestSetup = quote
    import AbstractAlgebra as AA
    using QuantumEnvelopingAlgebras
end
```

## 使い方

```jldoctest
julia> import AbstractAlgebra as AA

julia> using QuantumEnvelopingAlgebras

julia> R, q = AA.rational_function_field(AA.QQ, :q);

julia> U, (E, F, K, Ki) = uqsl2(R, q);

julia> E*F - F*E == (K - Ki)/(q - inv(q))
true

julia> K*Ki == one(U)
true

julia> iszero(commutator(casimir(U), E))
true
```

正規形は `F^a K^b E^c`（`a, c ≥ 0`, `b ∈ ℤ`）の有限和です。`Ki` は `K^(-1)` を
表します。`uqsl2(R, q)` は呼び出しごとに新しい親オブジェクトを作り、異なる親の元
どうしは演算できません。

`pbw_terms` は `(PBWMonomial, 係数)` の組を `(a+c, a, b, c)` の降順で返します。

```jldoctest
julia> R, q = AA.rational_function_field(AA.QQ, :q);

julia> U, (E, F, K, Ki) = uqsl2(R, q);

julia> pbw_coefficient(E*F, 0, 1, 0)
q//(q^2 - 1)
```

`uqsl2(R, q; max_terms=1000)` と指定すると、各演算段階の正規形に保持する PBW 項数
を制限できます。超過時は [`ComputationLimitError`](@ref) を投げ、不完全な元は返し
ません。

## Hopf 構造

```jldoctest
julia> R, q = AA.rational_function_field(AA.QQ, :q);

julia> U, (E, F, K, Ki) = uqsl2(R, q);

julia> T2 = tensor_power(U, 2);

julia> Delta = coproduct_map(U, T2);

julia> Delta(E) == tensor(T2, E, one(U)) + tensor(T2, K, E)
true

julia> antipode(E) == -Ki*E
true
```

## 量子平面

`XY = qYX` を満たす量子平面の正規形は `Y^a X^b`（`a, b ≥ 0`）です。

```jldoctest
julia> R, q = AA.rational_function_field(AA.QQ, :q);

julia> P, (X, Y) = quantum_plane(R, q);

julia> X*Y == q*Y*X
true
```

```@contents
Pages = ["index.md", "api.md"]
```
