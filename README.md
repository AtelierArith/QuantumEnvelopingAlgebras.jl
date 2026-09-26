# QuantumEnvelopingAlgebras.jl

`U_q(sl_2)` と量子平面を厳密な PBW 正規形で計算する Julia パッケージです。係数体には AbstractAlgebra.jl を使います。

```julia
import AbstractAlgebra as AA
using QuantumEnvelopingAlgebras

R, q = AA.rational_function_field(AA.QQ, :q)
U, (E, F, K, Ki) = uqsl2(R, q)

E*F - F*E == (K - Ki)/(q - inv(q)) # true
K*Ki == one(U)                      # true
commutator(casimir(U), E) == zero(U) # true

pbw_terms(E*F)                      # sorted (PBWMonomial, coefficient) pairs
pbw_coefficient(E*F, 0, 1, 0)      # 1/(q - q^(-1))
```

正規形は `F^a K^b E^c`（`a,c ≥ 0`, `b ∈ ℤ`）の有限和です。`Ki` は `K^(-1)` を表します。`uqsl2(R, q)` は呼び出しごとに新しい親オブジェクトを作り、異なる親の元どうしは演算できません。

現在の対象は `QQ(q)` と `QQ` 上の `q ≠ 0, ±1` です。非負整数べき、係数による除算、`s*K^b` の逆元と負べきを扱います。一般の非可換元の逆元は未対応です。

`uqsl2(R, q; max_terms=1000)` と指定すると、各演算段階の正規形に保持する PBW 項数を制限できます。超過時は `ComputationLimitError` を投げ、不完全な元は返しません。

## Hopf 構造

```julia
T2 = tensor_power(U, 2)
Delta = coproduct_map(U, T2)

Delta(E) == tensor(T2, E, one(U)) + tensor(T2, K, E) # true
counit(K) == one(R)                                   # true
antipode(E) == -Ki*E                                  # true
```

`tensor_power(U, n)` は毎回新しい親を作ります。`tensor(T2, ...)` のように演算先の親を明示してください。`coproduct_on_slot(T3, x, i)` はテンソル元 `x` の第 `i` スロットに余積を適用し、平坦化した `T3` の元を返します。

```julia
rho = highest_weight_representation(U, 3)
rho(E) * rho(F) == rho(E*F) # true
```

`rho(x)` は AbstractAlgebra の `4 × 4` 行列です。

```julia
phi = specialization_map(U, AA.QQ, AA.QQ(2))
phi(E*F) == phi(E)*phi(F) # true
```

特殊化は `QQ(q)` の各係数を約分済みの有理関数として評価します。分母が特殊化点で零になる元に `phi` を適用すると `SpecializationError` になります。`phi.target` は写像に固有の親オブジェクトです。

## 量子平面

第二の PBW 代数として、`XY=qYX` を満たす量子平面を実装しています。正規形は `Y^aX^b`（`a,b ≥ 0`）です。[この関係式を用いる文献](https://eprints.whiterose.ac.uk/id/eprint/126169/1/connqweyl.pdf)と同じ向きを採用しています。

```julia
P, (X, Y) = quantum_plane(R, q)
X*Y == q*Y*X            # true
basis_monomial(P, 2, 3) # Y^2*X^3
```

`q=1` は可換な場合として扱い、`q=0` は拒否します。量子平面の親も構築ごとに固有です。`quantum_plane(R, q; max_terms=1000)` で項数上限を指定できます。

```julia
import Pkg
Pkg.test()
```

簡単な所要時間と項数の計測は次のように実行します。

```sh
julia --project=. benchmark/benchmarks.jl
```

数学的規約と今後の範囲は [仕様書](docs/agents/QuantumEnvelopingAlgebras-spec.md) を参照してください。

## 開発

このパッケージの実装とドキュメントの整備には、OpenAI Codex と DeepSeek v4.1 を使用しました。

## License

このパッケージは [GNU General Public License v3.0](LICENSE)（GPL-3.0）の下で配布されます。AbstractAlgebra.jl は別のライセンス（個別ソースは BSD-2-Clause を含む）で配布されます。
