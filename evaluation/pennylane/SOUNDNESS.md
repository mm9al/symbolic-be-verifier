# PennyLane 優化的誤差不變量與檢查範圍

本文件對應 `bounded_pruning.py`、`local_summaries.py`，說明目前的數學論證、
實作所追蹤的對象，以及 exact algebraic checks 與 formal verification 的區別。
這些檢查沒有經過 proof assistant 核驗，不能將整個 rewriter 稱為 formally verified。

## 1. 誤差涵蓋所有輸出 ancilla branches 與任意 system input

令 `J_0 = |0^a> ⊗ I_S`，其中 `I_S` 涵蓋整個 system Hilbert space。
在第 t 個 gate 後，原始電路所對應的算子是

```text
V_t = G_t ... G_1 J_0 = sum_b |b> ⊗ B_(t,b).
```

`BranchState` 儲存此 rectangular operator 的全部輸出 branches，並非只儲存
某個 system state 的振幅，也不是只儲存最後要投影出的 `b=0` branch。
ancilla input 固定為零；目前沒有追蹤完整 square unitary 的所有 ancilla input columns。

考慮一次 pruning 的改動：

```text
D_t = sum_(b,P discarded) d_(b,P) |b> ⊗ P.
||D_t||_op <= sum_(b,P discarded) |d_(b,P)|.
```

每個 `|b> ⊗ P` 是從 system space 到 joint ancilla/system space 的 norm-1 算子。
程式以所有 `state.branches.items()` 中被刪除係數的 `abs_upper()` 累加並向上取界，
實作此不等式。即使當下的零 branch 完全未變，其他 branch 的成本仍會記帳。
後續 gates 可以把這些 branch 混回最後選取的 branch；不能只看當下的 top-left block。

更精確地說，若保留的 Acb balls 尚有半徑，不變量是：**存在**一個與所有保留係數球
相容的算子 `W_t`，使 `||V_t-W_t||_op <= eta_t`。不宣稱每個球的任意 midpoint
都自動滿足此界。成功 recenter 的係數則已換成記帳後的精確 dyadic midpoint。

所有支援的後續 gates 都是 joint-space unitary，所以

```text
||G_t(V_(t-1)-W_(t-1))||_op = ||V_(t-1)-W_(t-1)||_op.
```

區間運算包住對 `W_(t-1)` 的精確 gate action；新的 pruning／recentering 改動再依
三角不等式加進 eta。若一批操作會超過預算，就保留整批原值及原有區間。
最後投影 `(<0^a| ⊗ I_S)` 的 norm 是 1，因此原始 block 與計算 block 的差也被 eta 覆蓋。
還原的 global phase 模長為 1，數值求值的剩餘不確定性則由最後的係數球包住。

這足以支撐目前固定 ancilla input 的 block-encoding 驗證。
**沒有據此宣稱 `||U-U_tilde||_op <= eta` 對所有 ancilla input 同時成立**。
若要提供 full-unitary approximation certificate，需要另追蹤所有 input columns，
或另外建立完整 square operators 的誤差界。

## 2. Acb 收窄界住的是複數距離

設係數球包含

```text
z = (m_R + e_R) + i(m_I + e_I),
|e_R| <= r_R, |e_I| <= r_I.
```

換成 midpoint `m_R+i*m_I` 的誤差滿足

```text
|z-(m_R+i*m_I)| <= sqrt(r_R^2+r_I^2) <= r_R+r_I.
```

`ErrorLedger.recenter` 使用
`(ball.real.rad() + ball.imag.rad()).upper()`，再對所有 branches 的所有被收窄係數
向上累加，沒有只使用其中一個半徑。`ball.mid()` 是精確 dyadic 中心；不能另外轉成
binary64 中心卻不記錄轉換誤差。若未收窄，原始 balls 繼續保留其半徑。

最後 acceptance 使用包含保留區間半徑的 rigorous residual upper bound R，並要求

```text
R = 2^(n/2) * ||alpha*c_tilde - h||_2 upper bound,
alpha > 0,
R + |alpha|*eta <= epsilon_BE.
```

故原始 block 滿足 `||alpha*B-H||_op <= epsilon_BE`。
Hermiticity 的 coefficient tolerance 也包含 eta：
`max_P |Im(c_tilde_P)| + eta < delta`；此規則不代表精確 Hermiticity，
也不代表整個 anti-Hermitian operator 的 norm 小於同一個 delta。

## 3. Local summary 的 identity 與可信範圍

`certify_template` 對每個三量子位 template，以 `Q(sqrt(2), i)` 中的精確算術，
比較全部 8 個 computational-basis input columns。
role 0 是 basis bit 的最高位；roles 的順序是 control 0、control 1、target。
gates 按原始時間順序套用，即比較 `G_L ... G_1` 與 CCX／CCZ 的矩陣。

比較的是完整複數振幅，沒有消除、最佳化或忽略 global phase。
因此相差 `-1`、`i` 或 `exp(i*pi/4)` 的 template 不能通過。
8 個 basis columns 全部相等，依線性性可推出這個有限 template 的完整矩陣 identity，
包含與外部系統糾纏的輸入；它不是只比較某幾個數值樣本。

`_match` 要求 gate 順序、gate 名稱、arity、每個 operand 的 role mapping 一致，
並要求三個 roles 映射到三個不同 wires。傳回的 tuple 保留 role 順序，沒有排序 wires。
將 identity 嵌入其他 wires，相當於對兩邊同時 tensor identity 並作相同的 wire permutation，
故矩陣相等仍成立。`apply_summary` 再透過 `ancillas.index`／`systems.index` 對應
實際 branch bit 與 Pauli-string 位置，不以實際 wire number 大小推測位元順序。

需要區分三層依據：

- **數學論證**：上述 norm 不變量、矩陣相等的線性性與 wire permutation 的嵌入性質。
- **精確代數檢查**：Python/SymPy checker 對每個有限 template 全矩陣直接比較。
  JSON 中的 `certificates` 是此 checker 的檢查紀錄，並非獨立 proof kernel 可驗證的證明物件。
- **Regression tests**：檢查 matcher、更新 kernel、phase、wire mapping、dirty ancillas、
  區間半徑與誤差預算的實作。通過測試不能推出整個 rewriter 已形式化驗證。

目前仍信任 Python 實作、SymPy 精確算術、Arb/Acb 的 enclosure 行為、parser、
primitive gate semantics、matcher 與 summary execution 的正確性。
適合的描述是 **exact algebraically checked local identities with regression-tested
matching and execution**，而非 formally verified rewrite system。
