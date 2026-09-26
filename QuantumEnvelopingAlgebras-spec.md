# QuantumEnvelopingAlgebras.jl — 実装仕様書

- 仕様版: 0.1-draft
- 作成日: 2026-09-25
- 状態: 設計提案。実装状況は `README.md` と `src/` を参照。
- 仮称: `QuantumEnvelopingAlgebras.jl`。General registry の名称重複確認は未実施。
- 新規に作成するソースコードのライセンス方針: MIT。
- 最初の計算対象: Drinfeld–Jimbo 型の量子包絡代数 `U_q(sl_2)`。
- 文中の `[S01]` 等の出典は `REFERENCES.md` を参照。

## 1. 決定事項

AbstractAlgebra.jl を係数体と代数インターフェースの基盤に採用する。初版は任意の生成関係を解く CAS ではなく、既知の PBW 基底を持つ `U_q(sl_2)` の厳密演算パッケージとする。

本体の元は `F^a K^b E^c` の疎な線形結合として保存する。`a,c >= 0`、`b` は整数。掛け算のたびに正規形にする。一般的な式木を保持して後から `simplify` する設計は採らない。

非可換 Gröbner 基底の自動計算は実験・照合経路に限定し、通常の積・等号判定の前提にしない。AbstractAlgebra.jl に自由結合代数および非可換 Gröbner 基底計算が存在する点は確認済みである。ただし上限によって打ち切られた結果は完全な基底とは限らない。[S03, S04]

## 2. 既存基盤の評価

### 2.1 AbstractAlgebra.jl

採用するのは次の機能である。

- 有理数体、有理関数体と係数の厳密演算。
- `NCRing` / `NCRingElem`、parent と element の分離。
- 自由結合代数と、試験的な生成関係の表現。

独自実装するのは `U_q(sl_2)` の親オブジェクト、PBW 元、正規化された積、必要に応じた tensor/Hopf 構造である。AbstractAlgebra の商環 API だけでこれらがすべて自動的に実現する、とは仮定しない。[S03–S07]

### 2.2 OSCAR.jl

OSCAR は PBW 代数、イデアル、PBW 代数の商である GR 代数を提供する。汎用イデアル計算を目的にする場合には有力である。`pbw_algebra` の基底条件チェックも存在する。[S08, S09]

一方、このプロジェクトでは直接依存にしない。理由は、独自の MIT コアを小さく保つことと、最初の対象である `U_q(sl_2)` の正規形を直接実装できるためである。OSCAR の `.jl` コードを MIT 本体へ転用しない。OSCAR は GPL-3.0-or-later である。[S10]

`K^{-1}` は通常の非負指数の PBW 変数と同じではない。OSCAR 側でモデル化するなら、逆元用変数と逆元関係を入れた商などを別途検討する。本仕様は動作確認済みの OSCAR コンストラクタを提供するものではない。

### 2.3 QuaGroup

GAP の QuaGroup は量子包絡代数、Hopf 構造、最高ウェイト加群などの比較対象になる。ただしライセンスは GPL-2.0-or-later で、本体依存にはしない。[S11]

比較実験を将来設ける場合にも、規約の変換とライセンスの扱いを個別に確認する。別プロセス化や optional dependency 化だけで配布条件が自動的に解消するとは判断しない。

## 3. リリース範囲

### 3.1 v0.1: 代数としての正確な計算

必須機能:

- 明示的な係数体 `R` と変形パラメータ `q` を指定した構築。
- 生成元 `E,F,K,Ki`、`Ki = K^{-1}`。
- ゼロ、単位元、加減算、積、係数との演算、非負整数べき。
- 証明可能な形 `s*K^b` (`s != 0`) に対する逆元・負べき。
- PBW 正規形、係数取得、交換子、Casimir 元、決定的な表示。
- 親オブジェクト検証、正確な等号判定、copy/deepcopy/hash の契約。
- 独立な語の書き換えと行列表現による試験。

初版の保証対象は `QQ(q)` および `QQ` 上の `q0 != 0, +/-1` とする。コンテナ型は他の厳密な標数 0 の体も表現できる設計にするが、未試験の体への対応を宣言しない。

### 3.2 v0.2: Hopf 構造と表現

- `U` の有限テンソル冪、余積、余単位、antipode。
- 型 1 の有限次元最高ウェイト表現。
- 分母を検査する明示的な係数特殊化。

### 3.3 それ以降

量子平面・Weyl 代数などの第二の実装を加えてから共通 PBW カーネルを抽出する。ユーザー定義の生成関係、非可換 Gröbner 基底との橋渡し、高ランク量子群は別の設計段階とする。

初版の対象外:

- 任意の有限表示代数について常に停止する等号判定。
- 任意の非可換イデアルの Gröbner 基底計算を通常演算の一部として実行すること。
- `U_q(sl_n)` の一般実装、Lusztig integral form、divided powers、small quantum group。
- 普遍 R 行列を通常の有限疎多項式として格納すること。
- `q=1` や `q=-1` への機械的な代入、浮動小数による記号等号判定。
- `adjoint(E)=F` などの * 構造を暗黙に導入すること。

## 4. 数学的定義と規約

### 4.1 係数体

既定の使用例は `R = QQ(q)` とする。`q` は `U` の非可換生成元ではなく、係数体の元であり、すべての代数元と可換である。[S05]

```julia
import AbstractAlgebra as AA
R, q = AA.rational_function_field(AA.QQ, :q)
```

コンストラクタは次を検査する。

1. `R` がサポート対象の厳密な標数 0 の体である。
2. `q` を `R` に明示的に変換できる。
3. `q != 0`。
4. `q - inv(q) != 0`。

一般の環まで拡張する将来版では、3 と 4 の「非零」は不十分であり、両者の可逆性を検証する必要がある。

### 4.2 生成関係

`d = q - q^{-1}` とする。

\[
KK^{-1}=K^{-1}K=1,\qquad
KEK^{-1}=q^2E,\qquad
KFK^{-1}=q^{-2}F,
\]
\[
EF-FE=\frac{K-K^{-1}}{d}.
\]

この規約は [S12] Definition 1.1 と一致する。文献による `q` と `q^2` の違い、および余積の反転を混在させない。

### 4.3 正規形

\[
M(a,b,c)=F^aK^bE^c,\qquad a,c\in\mathbb N_0,\ b\in\mathbb Z.
\]

すべての元を有限和

\[
x=\sum_{a,b,c}\alpha_{a,b,c}M(a,b,c),\qquad\alpha_{a,b,c}\in R
\]

として保持する。正規形の一意性は付録 A の停止性・重複解消で裏づける。

`K` と `Ki` を別々の指数軸として本体に持たない。たとえば `K^3*Ki^5` は指数 `b=-2` にする。`Ki` は API 上の生成元だが、保存には独立変数を必要としない。

### 4.4 q 整数

\[
[n]_q=\frac{q^n-q^{-n}}{q-q^{-1}}
=\sum_{j=0}^{n-1}q^{n-1-2j}\quad(n>0),\qquad[0]_q=0.
\]

実装では有限和または等価な漸化式を使用できる。`[n]_q!` を割る計算は基本積に導入しない。

## 5. 内部データ構造

以下は設計スケッチであり、そのまま完全なパッケージになるコードではない。

```julia
struct PBWMonomial
    a::Int  # exponent of F, nonnegative
    b::Int  # exponent of K, may be negative
    c::Int  # exponent of E, nonnegative
end

mutable struct UqSl2Algebra{C,R<:AbstractAlgebra.Field} <: AbstractAlgebra.NCRing
    base::R
    q::C
    qinv::C
    dinv::C
end

struct UqSl2Elem{C,R<:AbstractAlgebra.Field} <: AbstractAlgebra.NCRingElem
    parent::UqSl2Algebra{C,R}
    terms::Dict{PBWMonomial,C}
end
```

`Ring` / `RingElem` ではなく `NCRing` / `NCRingElem` を用いる。AbstractAlgebra の型階層では前者は可換側の型である。[S03, S06]

`C<:FieldElem` と安易に制限しない。AbstractAlgebra は Julia 標準の `Rational` 等も係数として使用するためである。`C = elem_type(R)` とし、構築時に型と係数 parent を検証する。[S06]

`UqSl2Algebra` が mutable であるのは一意な parent identity を持たせるためであり、公開後に `q` や規約を書き換えてよいという意味ではない。親の構成情報は公開 API 上で不変とする。

### 5.1 不変条件

- すべてのキーで `a,c >= 0`。
- 係数の型と所属先は `R` と一致する。
- 零係数の項を保存しない。
- 同一単項式の係数は集約する。
- 空の辞書を零元とする。
- 辞書の参照を外部に公開しない。
- 外部から渡された可変係数と内部保存領域を共有しない。
- 通常演算は入力を変更しない。

### 5.2 親オブジェクト

初版の `uqsl2(R,q)` は新しい parent を作る。同じパラメータで構築しても、異なる parent の元は自動混合しない。

- `x+y`, `x*y`: parent 不一致なら `ParentMismatchError`。
- `x==y`: 異なる parent の代数元どうしなら `false`。
- `U(x)`: `parent(x) === U` の場合のみ同じ代数元として受理。
- 係数は `R` からの中央埋め込みで受理。
- 同型な代数間の移送は、将来の明示的な写像に限定する。

parent の自動キャッシュは v0.1 では設けない。v0.2 のテンソル親には明示的な親オブジェクトを API に渡し、親の相違で同じ式が加算不能にならないようにする。

### 5.3 copy、等号、ハッシュ

`deepcopy(x)` は係数と保存領域をコピーするが、parent はコピーせず元と同一に保つ。[S07]

代数元どうしの `isequal` は同一 parent と同じ正規化係数を要求する。ハッシュは parent identity と正規化項から作り、辞書の挿入順序に依存させない。`isequal(x,y)` なら必ず `hash(x)==hash(y)` とする。

`x == s` の係数比較を提供する場合は `x == U(s)` として実装する。一方、異種型間の `isequal(x,s)` は初版では `false` とし、parent 依存ハッシュとの整合性を壊さない。この違いは API 文書に明記する。

表示順は有限項を `(a+c,a,b,c)` の辞書式降順に整列する。この表示規則は Gröbner 基底の admissible order を主張するものではない。`a+c` を含む整数計算はオーバーフローを検査する。

## 6. 積の実装

### 6.1 基本方針

任意の積を文字列や一般 Expr の置換で処理しない。PBW 単項式に対する右乗算を実装し、分配法則で疎和に拡張する。

### 6.2 E と K の右乗算

\[
M(a,b,c)E^t=M(a,b,c+t),\qquad t\ge0,
\]
\[
M(a,b,c)K^t=q^{-2ct}M(a,b+t,c),\qquad t\in\mathbb Z.
\]

`K^t` は `abs(t)` 回のループにせず、指数更新と係数の積で処理する。

### 6.3 F の右乗算

生成関係から帰納的に

\[
E^cF=FE^c+\frac{[c]_q}{d}
(q^{1-c}K-q^{c-1}K^{-1})E^{c-1}
\]

を得る。従って `c>0` のとき

\[
\begin{aligned}
M(a,b,c)F={}&q^{-2b}M(a+1,b,c)\\
&+\frac{[c]_q q^{1-c}}{d}M(a,b+1,c-1)\\
&-\frac{[c]_q q^{c-1}}{d}M(a,b-1,c-1).
\end{aligned}
\]

`c=0` のときは第一項のみ。負の E 指数を一度構築してから零係数で消す実装は禁止する。

特に、第一項の `q^(-2b)` を後二項へ一律に掛けない。この係数は K と新しい F を交換するときだけ生じる。

### 6.4 単項式どうしの積

右側の単項式が `M(r,s,t)=F^r K^s E^t` なら、左側の単項式に対して

1. `right_mul_F` を `r` 回適用。
2. `right_mul_Kpow` を指数 `s` で一度適用。
3. `right_mul_Epow` を指数 `t` で一度適用。

する。各段階で同じキーを集約し、零を除去する。一般の元の積は項の組ごとの結果を集約する。

`E^c F^r` では k 回の交換子を使った項の E/F 指数はそれぞれ `c-k,r-k`、K 指数は `k,k-2,...,-k` となる。したがって単項式積の異なるキー数は、高々 `(m+1)(m+2)/2`、`m=min(c,r)` である。これは出力支持の上界であって、係数のビット複雑性や計算時間の上界ではない。

### 6.5 最適化の順序

まず上記漸化式の単純な実装を完成させる。次に q の整数べきと q 整数の再利用、`E^c F^r` の動的計画法、ビルダーの割り当て削減をベンチマークに基づき導入する。

グローバルで無制限の式キャッシュを作らない。キャッシュは parent や明示的 workspace に所属させ、上限とスレッド安全性を定める。

## 7. v0.1 公開 API

以下は新設する API の契約であって、既存パッケージにすでに存在する関数の案内ではない。

```julia
U, (E,F,K,Ki) = uqsl2(R, q)

parent(E) === U
base_ring(U) === R
quantum_parameter(U) == q

gens(U)  # (E,F,K,Ki), この順序
ngens(U) # 4。PBW 指数の軸数 3 とは異なる。

zero(U)
one(U)
U(2)
U(q + inv(q))

basis_monomial(U, a, b, c)  # F^a K^b E^c
pbw_coefficient(x, a, b, c) # 欠けている項は zero(R)
pbw_terms(x)                # ソート済みの独立したスナップショット
normal_form(x)              # x は常に正規形。値を変えない。
commutator(x, y)             # x*y-y*x
qinteger(q, n)              # n >= 0
casimir(U)
```

型インターフェースとして `parent`, `base_ring`, `elem_type`, `parent_type`, `base_ring_type`, `characteristic`, `is_exact_type` を確認・実装する。`zero`, `one`, `iszero`, `isone`, arithmetic, display, deepcopy, equality, hash を試験する。[S06, S07]

可換環の `gcd`, `factor`, 一般の分数体構築などがこの型で使える、と宣言しない。AbstractAlgebra の型階層に参加することと、全ての可換環アルゴリズムの入力になれることは異なる。

### 7.1 逆元・除算

`inv(x)` は `x=s*K^b, s!=0` と検証できる場合のみサポートし、結果は `inv(s)*K^(-b)` とする。

一般の非可換元を分母にした `/` は実装しない。`x/s` は `s` が係数体の非零元の場合に限り `inv(s)*x` とする。サポート外の逆元は `UnsupportedInverseError` とし、「調べていない」ことと「数学的に逆元が存在しない」ことを混同しない。

`x^n, n>=0` は繰り返し二乗法。負べきは上記のサポートされた逆元に限定する。`zero(U)^0=one(U)` とする。指数の符号反転を含め、`typemin(Int)` などのオーバーフローを検査する。

### 7.2 Casimir 元

本仕様では次の非スケール版を `casimir(U)` とする。

\[
C=FE+\frac{qK+q^{-1}K^{-1}}{d^2}.
\]

[S13] の Casimir は本仕様の C を `d^2` 倍した規約なので、比較時には変換する。`[C,E]=[C,F]=[C,K]=0` を厳密に試験する。

### 7.3 エラー

| 状況 | 挙動 |
|---|---|
| 非厳密係数体、非対応の体 | `ArgumentError` |
| q=0、q-inv(q)=0 | `DomainError` |
| 別 parent の元を演算 | `ParentMismatchError` |
| a<0 または c<0 | `DomainError` |
| Int 演算の範囲超過 | `OverflowError` |
| 未対応の逆元 | `UnsupportedInverseError` |
| 特殊化で分母が零 | `SpecializationError`（v0.2） |
| 明示的な計算予算を超過 | `ComputationLimitError` |

独自例外はパッケージで定義する。計算上限を超えた場合に項を黙って切り捨てたり、不完全な正規形を正規化済みの元として返したりしない。

## 8. v0.2: Hopf 代数仕様

余積の規約は [S13] と一致させる。

\[
\Delta(K)=K\otimes K,\quad
\Delta(E)=E\otimes1+K\otimes E,\quad
\Delta(F)=F\otimes K^{-1}+1\otimes F.
\]
\[
\epsilon(K)=1,\quad\epsilon(E)=\epsilon(F)=0,
\]
\[
S(K)=K^{-1},\quad S(E)=-K^{-1}E,\quad S(F)=-FK.
\]

`Delta` は代数準同型として、`S` は **反**準同型として拡張する。特に `S(x*y)=S(y)*S(x)` である。

### 8.1 テンソル積

`U` の有限テンソル冪を独立した `NCRing` として実装する。係数は同じ `R`、基底キーは `NTuple{N,PBWMonomial}` とする。積は成分ごとの積であり、異なるスロットの元は可換、同じスロットでは元の非可換積を用いる。

```julia
T2 = tensor_power(U, 2)
Delta = coproduct_map(U, T2)
Delta(E) == tensor(T2, E, one(U)) + tensor(T2, K, E)
```

これは提案 API である。`tensor(T2,...)` の明示的 parent により、同じ数学的対象の親オブジェクトの不一致を避ける。テンソルの括弧は内部的に平坦化し、余結合則の両辺を同じ `U^tensor3` で比較する。

必要な試験:

- `Delta(x*y)==Delta(x)*Delta(y)`。
- `(Delta tensor id)Delta(x)==(id tensor Delta)Delta(x)`。
- 左右の counit 恒等式。
- `m(S tensor id)Delta(x)==epsilon(x)*1` と右側の antipode 恒等式。
- 生成関係が `Delta` によって零へ移ること。

## 9. 行列表現と特殊化

### 9.1 最高ウェイト表現

標準試験では最高ウェイト n、基底 `v_0,...,v_n` を使う。

\[
Kv_j=q^{n-2j}v_j,
\]
\[
Fv_j=v_{j+1}\ (j<n),\qquad Fv_n=0,
\]
\[
Ev_j=[j]_q[n-j+1]_qv_{j-1}\ (j>0),\qquad Ev_0=0.
\]

これは `v_j=F^j v_0` の正規化であり、divided-power 基底ではない。[S12] の表現公式は基底スケーリングによってこの形にできる。行列の第 j 列を基底ベクトルの像とする。

`rho(F^a K^b E^c)=rho(F)^a rho(K)^b rho(E)^c` とする。順序を逆にしない。n=0,...,6 で生成関係と `rho(x*y)==rho(x)*rho(y)` を試験する。ただし有限個の行列表現で一致したことだけを、代数の等号判定の証明とはしない。

### 9.2 特殊化

`QQ(q)` の元の `q=q0` への代入は、**体全体に定義された準同型ではない**。たとえば `1/(q-q0)` には代入できない。

v0.2 の特殊化は、各元の各係数を約分済みの有理関数として評価し、すべての分母が非零である場合だけ定義する。ターゲットの `q0` 自体が生成関係の条件を満たすことも検査する。

`q=1` の古典極限は、適切な integral form や再パラメータ化を別に設計する必要があるため、通常の特殊化 API では拒否する。1 の冪根にパラメータを置き換えても、自動的に small quantum group になるわけではない。商関係や integral form の選択は別仕様とする。

## 10. テストおよび受入条件

### 10.1 生成関係

`K*Ki=Ki*K=1`、K と E/F の交換関係、E/F の交換子を `QQ(q)` 上で厳密に試験する。`K^(-m)` と `Ki^m` の一致も検査する。

### 10.2 独立な正規化器

`test/reference_word_rewriter.jl` に、本体の積を呼ばない語の書き換え器を置く。規則は付録 A の 7 本とする。記号係数で 8 個の重複を確認し、短い語について本体の PBW 積と照合する。

語の書き換え器を、`normal_form` と同じ本体関数の別名にしてはならない。

### 10.3 性質試験

固定乱数 seed を使い、結合則、左右の分配則、単位・零、係数の中央性、正規化の冪等性を試験する。負の K 指数を必ず含める。

異なる構築経路で同じ元を生成し、等号・ハッシュ・表示が一致することを確認する。親の相違、入力係数の alias、deepcopy での parent 保持も試験する。

### 10.4 行列表現・Casimir

9.1 の表現を本体の積と独立に構成し、生成関係、積の準同型性、Casimir の中心性を確認する。

### 10.5 受入条件

- `QQ(q)` 上の厳密試験に成功する。
- `QQ` 上の q=2,3/2 の試験に成功する。
- 文書中の v0.1 サンプルが CI で実行される。
- 上限超過・未対応値・parent 不一致が明示的なエラーになる。
- v0.1 の通常使用に Oscar/Singular/GAP のインストールが不要。
- 倍精度近似や行列表現だけに依存する等号判定を行わない。

## 11. リポジトリ構成案

```text
QuantumEnvelopingAlgebras.jl/
  Project.toml
  LICENSE
  THIRD_PARTY_NOTICES.md
  README.md
  src/
    QuantumEnvelopingAlgebras.jl
    parents.jl
    elements.jl
    coefficients.jl
    arithmetic.jl
    multiplication.jl
    printing.jl
    casimir.jl
  test/
    runtests.jl
    test_relations.jl
    test_parents.jl
    test_invariants.jl
    test_properties.jl
    test_representations.jl
    reference_word_rewriter.jl
  benchmark/
    Project.toml
    benchmarks.jl
  docs/
    Project.toml
    make.jl
    src/
      index.md
      conventions.md
      implementation.md
      limitations.md
```

v0.2 で `tensor.jl`, `hopf.jl`, `representations.jl`, `specialization.jl` を加える。

## 12. 依存関係と互換性

初期の直接 runtime 依存は AbstractAlgebra.jl のみとし、テストの `Test`, `Random` 等は適切に別宣言する。ベンチマークと文書の依存は別環境に置く。

```toml
[deps]
AbstractAlgebra = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"

[compat]
AbstractAlgebra = "0.50"
julia = "1.10"
```

これは互換性の初期提案であり、Julia で試験済みという主張ではない。参照した AbstractAlgebra v0.50.2 の Project.toml は Julia 1.10 を下限としている。[S02]

最初の CI は 1.10 とその時点の安定版を対象にする。パッケージ作成時に実際の manifest を生成し、採用版で関数名とディスパッチを検証する。未試験の 0.x マイナー版まで互換性を広げない。

PrecompileTools を導入する場合は直接依存として宣言する。まず正しさと初回実行時間を測り、小さな `QQ(q)` の構築と積だけを workload にする。大規模 Gröbner 基底、ファイル I/O、外部 CAS 起動を precompile workload に入れない。

## 13. 性能目標と安全な失敗

初版では根拠のないミリ秒目標を設定しない。`E^n F^n`, `(E+F)^n`、疎元どうしの積、負の K べき、QQ(q) 係数の分子分母成長について、時間・allocation・項数を記録する。

ボトルネックを非可換正規化と係数の有理関数演算に分けて測定する。浮動小数化による高速化は厳密演算の置き換えとして採用しない。

任意コード評価を伴う文字列入力は初版に導入しない。`eval` を用いる読み込み、任意 Julia オブジェクトの安全性を保証しないシリアライズは公開交換形式にしない。

## 14. ライセンス方針

新規実装を MIT とし、依存物を MIT へ再ライセンスしたと誤解させない。AbstractAlgebra の個別 `.jl` ソースは BSD-2-Clause であり、「AbstractAlgebra は MIT」と記載しない。[S01]

Julia 実行環境や依存ライブラリを同梱する配布物は、ソースパッケージ単体とは別にライセンスを確認する。AbstractAlgebra の LICENSE 自体も Julia 同梱ライブラリの扱いと個別ソースの扱いを分けている。[S01]

詳細は `LICENSE_NOTES.md` を参照。これは設計上の整理であり、個別の配布形態への法的判断を確定する文書ではない。

## 15. 高ランク化の境界

`U_q(sl_3)` 以降は、単純生成元 E_i/F_i を任意順に並べた指数辞書だけでは通常の PBW 基底の設計は完成しない。正根に対応した root vector、順序、量子 Serre 関係を仕様化する必要がある。

従って v0.1 を「n 個の変数へ拡張するだけで全量子群に対応」と説明しない。第二の代数を実装した段階で共通化できた演算だけを抽出し、一般化は段階的に行う。

## 付録 A. 語の書き換えによる正規形の確認

`L=K^{-1}` と置き、自由語上の規則を次の 7 本とする。

\[
KF\to q^{-2}FK,\quad LF\to q^2FL,\quad
EK\to q^{-2}KE,\quad EL\to q^2LE,
\]
\[
EF\to FE+d^{-1}K-d^{-1}L,\quad
KL\to1,\quad LK\to1.
\]

語長を優先し、同じ長さでは `E>K>L>F` の辞書順を使う。各規則の右辺の語は左辺より小さいので書き換えは停止する。

左辺は相異なる長さ 2 の語なので包含競合はなく、重複競合は以下の 8 個で尽くされる。左右いずれの位置から始めても次の同じ結果に還元される。

| 語 | 共通の還元結果 |
|---|---|
| EKF | `q^(-4)*F*K*E + q^(-2)*(K^2-1)/d` |
| ELF | `q^4*F*L*E + q^2*(1-L^2)/d` |
| EKL | `E` |
| ELK | `E` |
| KLF | `F` |
| LKF | `F` |
| KLK | `K` |
| LKL | `L` |

したがって既約語は `F^a K^b E^c` (`b>=0`) または `F^a L^r E^c` (`r>0`) に一意に還元され、整数 K 指数による本仕様の PBW 表現と一致する。

この確認は本仕様の規則についてのものであり、任意に追加した関係式について合流性を保証するものではない。

## 付録 B. この仕様作成時に実施した検算

`validation/reference_check.py` は Julia 本体とは独立の Python 参照計算である。SymPy と Fraction を用い、次を実行した。

- 8 個の重複を有理関数 q で記号的に検算。
- q=2 と q=3/2 のそれぞれで、長さ 0～6 の全 5,461 語について語の書き換えと PBW 漸化式を照合。
- 各 q で 100 組の結合則と分配則を検算。
- 各 q で 4 生成元との Casimir の交換を検算。
- 次元 1～7 の行列表現で生成関係を検算し、各 q で 35 組の積を照合。

結果は `validation/RESULTS.json` に保存した。これらは数学的規約・漸化式の検算であり、Julia の API、性能、実装のテストに代わるものではない。本環境では Julia 実行ファイルを確認できておらず、同梱 Julia サンプルは未実行である。


---

# 出典

確認日: 2026-09-25。設計提案と既存ソフトウェアの事実を区別するための一次資料一覧。

## 基盤・ライセンス

- [S01] AbstractAlgebra v0.50.2 LICENSE.md。個別 .jl ソースの BSD-2-Clause と、Julia 同梱ライブラリを含む配布物についての説明。
  https://raw.githubusercontent.com/Nemocas/AbstractAlgebra.jl/v0.50.2/LICENSE.md
- [S02] AbstractAlgebra v0.50.2 Project.toml。直接依存と Julia 互換性。
  https://raw.githubusercontent.com/Nemocas/AbstractAlgebra.jl/v0.50.2/Project.toml
- [S03] AbstractAlgebra: Free associative algebras。自由代数、NCRing 階層、語の取得、Gröbner 基底。
  https://nemocas.github.io/AbstractAlgebra.jl/stable/free_associative_algebra/
- [S04] 同資料の Gröbner bases 節。reduction_bound は追加された非零基底要素の数に関する上限であり、次数上限や wall-clock timeout ではない。不完全な基底への還元で非零でも、商で非零とは断定できない。
  https://nemocas.github.io/AbstractAlgebra.jl/stable/free_associative_algebra/#Groebner-bases
- [S05] AbstractAlgebra: Rational function fields。
  https://nemocas.github.io/AbstractAlgebra.jl/stable/function_field/
- [S06] AbstractAlgebra: AbstractTypes.jl および型インターフェース。
  https://raw.githubusercontent.com/Nemocas/AbstractAlgebra.jl/master/src/AbstractTypes.jl
  https://nemocas.github.io/AbstractAlgebra.jl/stable/types/
- [S07] AbstractAlgebra: Ring Interface。parent、係数変換、deepcopy、hash など。非可換型に可換専用操作まで実装するという意味ではない。
  https://nemocas.github.io/AbstractAlgebra.jl/stable/ring_interface/
- [S08] OSCAR: Creating PBW-Algebras。
  https://docs.oscar-system.org/stable/NoncommutativeAlgebra/PBWAlgebras/creation/
- [S09] OSCAR: GR-Algebras: Quotients of PBW-Algebras。
  https://docs.oscar-system.org/stable/NoncommutativeAlgebra/PBWAlgebras/quotients/
- [S10] Oscar.jl LICENSE.md。GPL-3.0-or-later。
  https://raw.githubusercontent.com/oscar-system/Oscar.jl/master/LICENSE.md
- [S11] QuaGroup 公式パッケージページおよび manual。GPL-2.0-or-later、量子包絡代数、Hopf 構造、加群。
  https://gap-packages.github.io/quagroup/
  https://gap-packages.github.io/quagroup/doc/chap3_mj.html
- [S14] Open Source Initiative: MIT License。
  https://opensource.org/license/mit

## 数学的規約

- [S12] T. Ito, P. Terwilliger and C.-w. Weng, “The quantum algebra U_q(sl_2) and its equitable presentation”, arXiv:math/0507477。Definition 1.1 の生成関係。Lemma 4.1 の表現は本仕様と基底スケーリングが異なる。
  https://arxiv.org/html/math/0507477
- [S13] Q. Labriet and L. Poulain d’Andecy, “Little q-Jacobi polynomials and symmetry breaking operators for U_q(sl_2)”, arXiv:2506.23848v1, Section 2.1。採用した coproduct と antipode が明示されている。Casimir のスケールは本仕様と異なる。
  https://arxiv.org/html/2506.23848v1

## 出典を使う範囲

本仕様の PBW 右乗算漸化式、重複競合の表、実装アーキテクチャ、試験計画は上記の生成関係から独立に構成した。参照したソフトウェアの実装コードを転載・翻訳したものではない。

stable 文書や master の内容は将来変わりうる。パッケージ実装時には採用タグで API とライセンスを再確認し、実行したテスト環境を記録する。
