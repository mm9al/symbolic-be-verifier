# PennyLane Ising block-encoding evaluation (TASK 1)

本資料夾只跑 TASK 1。輸入來自 PennyLane `qml.PrepSelPrep`，與 RQ1 使用同一個
periodic transverse-field Ising Hamiltonian：

```text
H_n = sum_{i=0}^{n-1} Z_i Z_{(i+1) mod n} + sum_{i=0}^{n-1} X_i
```

J=g=+1。n=2 依既有 benchmark 慣例計入兩次週期邊，因此 H_2=2 Z_0 Z_1+X_0+X_1。
PennyLane 會合併這兩項；`n_terms` 記錄實際 LCU 項數 3，`target_lcu_terms` 記錄原始
項數 4。名義 alpha=2n，只作來源 metadata；接受時的 alpha 由 verifier 獨立推得。

## 執行

在 `symbolic-be-verifier` 目錄使用現有的兩個環境：

```bash
# 生成：PennyLane 環境。先用小範圍，另存輸入以保留既有實驗。
MPLCONFIGDIR=/private/tmp/research-pennylane-mpl \
XDG_CACHE_HOME=/private/tmp/research-pennylane-cache \
PYTHONDONTWRITEBYTECODE=1 ../cobble/.venv/bin/python \
  evaluation/pennylane/generate_benchmarks.py --sizes 2,3 \
  --out-dir /private/tmp/pennylane-ising-inputs

# TASK 1：verifier 環境。epsilon 直接是 epsilon_BE，不分配 TASK 2 預算。
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  evaluation/pennylane/run_block_encoding.py \
  --manifest /private/tmp/pennylane-ising-inputs/manifest.csv \
  --output-dir /private/tmp/pennylane-ising-task1 \
  --block-encoding-epsilon 1e-8 --hermiticity-delta 1e-12 \
  --precision 100 --timeout-sec 1800
```

生成預設尺寸與既有 RQ1 generator 一致：
`2,3,4,5,6,7,8,9,10,16,17,32,33,64,65,128,256`。
省略 `--sizes` 會生成完整尺寸；亦可用 `--sizes 2..8` 或 `--sizes 2,4,8,16`。
原生 gate 展開和具體 Pauli symbolic execution 的成本會隨尺寸增加。
runner 預設在第一個 timeout 後跳過更大的尺寸；`--continue-after-timeout` 可繼續。
語義 REJECT 會照常記錄並執行下一列。

只重跑已生成的子集合或寫到其他位置：

```bash
.venv/bin/python evaluation/pennylane/run_block_encoding.py \
  --manifest benchmarks/pennylane/delta_calibration_manifest.csv \
  --sizes 2,3 --output-dir /private/tmp/pennylane-ising-rerun
```

已保存的輸入索引分成三份：

- `delta_calibration_manifest.csv`：n=2,3,4 的固定校準子集合。
- `all_ising_manifest.csv`：目前已生成的 n=2..10,16,17,32,33,64,65,128；尚未包含 n=256。
- `manifest.csv`：最近一次 generator 的子集合，目前為 n=17,32,33,64,65。

重跑已保存的實驗時明確指定前兩份 manifest，避免受到最近一次生成子集合影響。
`--out-dir DIR` 可改生成位置；runner 使用 `--manifest DIR/manifest.csv`。
manifest 的 input paths 相對於 manifest 所在目錄，可搬移整個資料夾。
如需獨立建立生成環境，使用 Python 3.12 安裝本資料夾的 `requirements-generator.txt`；
生成程式不依賴 `cobble` 或 `external_framework` 的 Python 模組。

## 可選的精確局部摘要

`--local-summaries` 在 PennyLane evaluation 中啟用精確 CCX／CCZ 子電路摘要。
預設不啟用；原本 core verifier、自建 benchmark generator、RQ1／RQ2 evaluation
都沒有變更。不加此選項時沿用原本逐閘路徑及 CSV 欄位。

```bash
.venv/bin/python evaluation/pennylane/run_block_encoding.py \
  --manifest benchmarks/pennylane/all_ising_manifest.csv \
  --sizes 2..10,16 \
  --local-summaries \
  --block-encoding-epsilon 1e-8 --hermiticity-delta 1e-12 \
  --precision 100 --timeout-sec 1800 \
  --output-dir /private/tmp/pennylane-ising-local-summaries
```

啟用後預設輸出到 `evaluation/pennylane/results_local_summaries/`，因此不會覆寫
原本的 `results/` 或 `results_delta/`。`--output-dir` 仍可明確指定其他位置。
timeout 行為沿用原本設定：第一個 timeout 後跳過後續案例。

目前辨識 PennyLane 的兩種 15-gate 分解：Toffoli（CCX）與 CCZ。
角色可映射到任意 wire indices；兩個控制都必須是已宣告的 ancillas，target 可為
ancilla 或 system qubit。摘要對所有輸入有效，包括 target ancilla 已是 1 的情況。
不符合完整 pattern 或不支援的 layout 會保留原始 gates 執行。

每個首次使用的 template 都以 Q(sqrt(2),i) 的精確振幅，計算三個 qubit 的全部
8 個 computational-basis input columns，檢查與 CCX／CCZ 完全相等，包含絕對相位。
沒有數值 tolerance；精確代數檢查失敗時不使用摘要。每個 occurrence 都記錄原始 gate 範圍、
wire mapping、template hash 與檢查結果。
這是可信 Python／SymPy checker 的 exact matrix identity check；整個 matcher／rewriter
尚未經 proof assistant 驗證，regression tests 也不能充當形式化證明。

只有角度精確等於 0 的合法單量子位旋轉會移除。非零小數 RX／RY／RZ 均保留原值，
不改成理想 pi 分數、不截斷小係數，也不增加誤差預算；摘要的替換誤差為 0。
原始 QASM 與全域相位 metadata 不變，alpha、epsilon_BE、delta 判定沿用原規則。
不同傳播順序可能產生不同寬度的 Arb 區間，測試會檢查所有 branch 係數區間相容。
摘要等價性的依據是精確 template identity 檢查；區間相容是額外的回歸檢查。
global phase、wire ordering 與 checker 的可信範圍見 [SOUNDNESS.md](SOUNDNESS.md)。

JSON 的 `local_summaries` 保存 certificates 與 occurrences；CSV 額外記錄
`local_summary_count`、`summarized_source_gates`、`executed_steps`、
`removed_identity_gates`、`summary_compilation_sec`。`gate_count` 仍是原始 QASM 的 gate 數。
runtime 包含 template 證明、辨識和摘要執行。小數旋轉造成的 Pauli 項數膨脹仍可能存在，
較大尺寸的加速效果需逐例量測。

下表直接取自 `results_local_summaries/baseline/` 與 `results_local_summaries/` 中
目前保存的 CSV／JSON。兩者是不同批次的單次量測，使用相同 QASM、100 位、
epsilon_BE=1e-8、delta=1e-12、tracemalloc，每例皆使用新的 worker process。
未做 warmup／repetitions，因此加速比只是這兩次量測的比值。

<!-- BEGIN local-summary-comparison -->
| n | 原本逐閘 (s) | 局部摘要 (s) | 加速 | 原始 gates → 執行 steps |
|---|---:|---:|---:|---:|
| 2 | 1.279 | 1.258 | 1.0x | 53 → 20 |
| 3 | 5.837 | 0.242 | 24.1x | 427 → 177 |
| 4 | 35.785 | 3.077 | 11.6x | 978 → 482 |
<!-- END local-summary-comparison -->

上述 cases 的接受結果相同，具體 block 的係數區間皆相容。
局部摘要版另完成 n=5，7.784 秒、1462 gates → 592 steps、ACCEPT；
配對 baseline 資料夾只包含 n=2,3,4，因此不把 n=5 放入配對比較表。
精確數值見 [comparison.csv](results_local_summaries/comparison.csv)，
CSV／JSON、配對輸入／設定及係數區間的重檢結果見
[comparison_checks.json](comparison_checks.json)。早期回歸檢查紀錄另保留於
[verification_checks.json](results_local_summaries/verification_checks.json)。

測試：

```bash
.venv/bin/python -m pytest tests/test_pennylane_local_summaries.py -q
```

## 可選的稀疏傳播與誤差記帳

局部摘要消除了部分 Toffoli／CCZ 展開，但沒有消除小數旋轉產生的小係數。
既有 `results_local_summaries/` 的 n=5、6、7 分別保留 768、4096、16384 個
Pauli 項，而摘要後操作數只從 592 增加至 846、1182。每個後續 gate 都要處理這些項，
是 7.8s → 107.8s → 580.1s 的主要成本來源。n=6、7 的最後表示已達到完整的 4^n 項。
原本 delta 僅用於最後判定，不會在傳播時刪除任何係數。

新增模式只在 PennyLane evaluation 啟用，不更動 core 或 RQ1／RQ2：

```bash
.venv/bin/python evaluation/pennylane/run_block_encoding.py \
  --manifest benchmarks/pennylane/all_ising_manifest.csv \
  --sizes 2..10,16 --bounded-pruning \
  --pruning-threshold 1e-24 --pruning-budget 1e-14 \
  --block-encoding-epsilon 1e-8 --hermiticity-delta 1e-12 \
  --precision 100 --timeout-sec 1800 \
  --output-dir /private/tmp/pennylane-ising-bounded-pruning
```

`--bounded-pruning` 包含精確局部摘要。預設輸出到
`evaluation/pennylane/results_bounded_pruning/`，與原始結果、delta 校準及
`results_local_summaries/` 分開。原始 QASM、輸入 hash、角度、全域相位都保留。
第一個 timeout 後跳過後續案例。`--pruning-budget 0` 停用此模式的截斷與區間收窄，
仍使用局部摘要。

此模式對**所有 ancilla branches** 執行兩個可記帳的操作：

1. 只有係數的嚴格 `abs_upper` 不大於 `--pruning-threshold` 時才考慮刪除；
   每批成本是所有被刪除 `(branch, Pauli)` 係數的絕對值上界之和。
   區間含零並不是刪除的充分條件。
2. 為避免長電路反覆放大互不相關的 Arb/Acb boxes，將保留的係數球換成精確的
   dyadic midpoint 前，先把其 real radius + imaginary radius 全部加入誤差帳。
   這些數值不確定性透過 operator norm 傳播，不會無記帳地忽略。

兩者共用 `--pruning-budget`，預設固定為 1e-14，定義在未乘 alpha 的 block B。
任何一批會超支時，保留該批原值；不自動放寬門檻。門檻、預算與 epsilon、delta
是不同參數。可依應用設定，但每次接受都要重新檢查以下原始電路誤差界。

理由是 branch column `V = sum_a |a> B_a` 的每個 Pauli 項，其基底算子 norm 都是 1。
單次改動的 operator norm 因此不大於上述成本之和。後續完整 unitary 保持 operator norm，
最後 ancilla projection 不增加它；全域相位也不增加它。

此處 V 的 ancilla input 固定為零，但保留所有輸出 ancilla branches 與任意 system input。
它是完整 branch column 的界，不是完整 square unitary 對所有 ancilla inputs 的界。
區間仍有半徑時的不變量，以及 Acb 複數中心距離的推導，見 [SOUNDNESS.md](SOUNDNESS.md)。

累計上界 eta 給出：

```text
||B_original - B_tilde|| <= eta

alpha > 0
AND 2^(n/2) * ||alpha*c_tilde - h||_2 + alpha*eta <= epsilon_BE
AND (independent exact Hermiticity proof OR max_P |Im(c_tilde_P)| + eta < delta)
```

即使某個 Pauli 項已被刪除，其係數誤差仍被 eta 覆蓋。不能用稀疏結果中的零虛部
推論原始 block 精確 Hermitian；只有獨立結構證明，或在扣除 eta 後仍可排除零的
虛部證據，才提供相應的精確結論。alpha 是候選值，最後接受依據包含 `alpha*eta`。
`core_proportionality_pass` 僅針對稀疏結果，`proportionality_pass` 才包含原始電路的誤差帳。

每例 JSON 的 `bounded_pruning` 保存逐批成本、總上界、刪除／收窄數量、超支跳過數與
peak/final 項數。`actual_block_coefficients` 在此模式表示 `B_tilde`，其範圍另記於
`coefficient_scope`；原始 block 的關係由 eta 證書給出。CSV 的
`truncation_error_upper` 是 pruning 與 charged recentering 的總 eta；
`roundoff_error_upper` 另列收窄區間的部分，`operator_error_upper` 已包含 `alpha*eta`。
未啟用此模式時，既有 CSV 欄位與計算路徑保持原樣。

```bash
.venv/bin/python -m pytest tests/test_pennylane_bounded_pruning.py -q
```

測試涵蓋全 branches 誤差界、非交換 unitary suffix、區間半徑記帳、預算用盡保留原值、
alpha 誤差縮放、Hermiticity 嚴格邊界、原始 QASM hash 與不同模式的輸出隔離。

下表使用 `results_local_summaries/` 與 `results_bounded_pruning/` 目前保存的
CSV／JSON，設定為 100 位、epsilon_BE=1e-8、delta=1e-12、threshold=1e-24、budget=1e-14。
兩欄來自不同批次的單次量測；配對案例的輸入 hash、目標與精度／門檻皆相同。
未做 repetitions 或統計平均。TIMEOUT／SKIPPED_AFTER_TIMEOUT 沒有可比較的完成時間。

<!-- BEGIN bounded-pruning-comparison -->
| n | 既有局部摘要 (s) | 新模式 (s) | 最後 Pauli 項數：既有 → 新模式 |
|---|---:|---:|---:|
| 5 | 7.784 | 1.239 | 768 → 70 |
| 6 | 107.801 | 3.026 | 4096 → 121 |
| 7 | 580.148 | 8.451 | 16384 → 174 |
| 8 | TIMEOUT | 28.772 | — → 239 |
| 9 | SKIPPED_AFTER_TIMEOUT | 29.896 | — → 269 |
| 10 | SKIPPED_AFTER_TIMEOUT | 35.272 | — → 361 |
| 16 | SKIPPED_AFTER_TIMEOUT | 385.061 | — → 1019 |
<!-- END bounded-pruning-comparison -->

本次 n=2..10、16 全部 ACCEPT。n=5..10 的總 eta 約 2e-28 至 3e-27，
n=16 的總 eta 約 1.59e-26、operator error upper 約 2.00e-11、最大虛部上界約 2.57e-16。
所有接受結果都包含誤差帳，epsilon_BE、delta 與 precision 未變。
原始數值、輸入與係數區間相容性檢查見
[comparison.csv](results_bounded_pruning/comparison.csv)。區間相容是回歸檢查，
誤差證書的依據仍是逐步的 norm 上界。
[validation_checks.json](results_bounded_pruning/validation_checks.json) 記錄原有 core／RQ1／RQ2
來源 hash 不變，以及每個測量電路的 unitary operands 檢查。
`source_snapshot/` 中的兩個模式檔案均與目前這批量測 metadata 的來源 hash 相符；
逐檔驗證結果見 `source_snapshot_checks`。歷史 `sources_sha256` 保留實際量測版本，
不會改寫成目前工作目錄的版本。
稀疏化仍保留顯著的小數誤差項，不能據此推論所有更大尺寸都會快速完成。

另一批完整 sweep 保存在 `results_bounded_pruning_all_ising/`：n=2..10,16,17 ACCEPT，
n=32 在 1800 秒 TIMEOUT，n=33,64,65,128 因 timeout 跳過。
其時間不混入上述固定比較批次。

重新產生比較表、README 表格及檢查紀錄，不重跑電路：

```bash
.venv/bin/python evaluation/pennylane/refresh_comparisons.py
.venv/bin/python evaluation/pennylane/refresh_comparisons.py --check
```

工具會先檢查 CSV 與案例 JSON 是否屬於同一批輸出、配對輸入／設定是否一致、
係數區間是否相容，以及保存的 source snapshots 是否符合量測 hash；不一致時停止。
`comparison_checks.json` 記錄產生比較表所使用的 CSV、JSON 與 run metadata hashes。

## 原生輸入與接受條件

直接分解同一個 PennyLane `PrepSelPrep` 物件，保留原生 PREPARE／SELECT／unprepare
segments，再匯出為 verifier 支援的 H/X/Y/Z/S/T/RX/RY/RZ/CX/CY/CZ basis。
沒有重建 SELECT、沒有 snapping decimal pi、沒有新增 workspace。
PennyLane 的無條件 GlobalPhase 逐項保存於 `structure.json`；其精確 decimal 總和
gamma 以 `exp(-i*gamma)` 乘回實際 block。受控 relative phases 留在分解後的 gates 中。

TASK 1 使用：

```text
(Exact Hermiticity OR max_P |Im(c_P)| < delta)
AND alpha > 0
AND ||alpha*c - h||_2 <= epsilon_BE / 2^(n/2)
```

`c`、`h` 是實際 block 與獨立 Ising target 的 Pauli 係數。
使用既有 verifier 的 symbolic execution 與 proportionality fit；選取具體 decimal
alpha 後，再以 Arb 重新計算殘差及精確 epsilon_BE 門檻。
精度至少 max(100, requested_precision, 2n) decimal digits。

`--hermiticity-delta` 預設 `1e-12`。此 delta 定義在**最後還原全域相位後、尚未乘上 alpha 的
block B 的 Pauli 係數**，取所有係數虛部絕對值的最大值。
比較使用 Arb 的嚴格上界，要求上界 **< delta**；相等不接受，也不以 midpoint 取代區間。
proportionality 的 epsilon_BE 檢查仍獨立決定是否通過。

目前校準資料為 n=2,3,4，n=3 的最大虛部上界約 `4.83e-16`，n=4 約 `5.58e-16`。
先用實測最大值的 1000 倍，再向上取十的冪，選定固定預設值 `1e-12`。
這是小型案例的初始設定；後續稀疏模式已量測 n=2..10、16，結果見上表，
沒有依每一列的結果提高 delta。
數值與輸入 provenance 見 [hermiticity_delta_calibration.json](hermiticity_delta_calibration.json)。
門檻不會依每一列的結果自動提高。可用 `--hermiticity-delta 0` 恢復原本精確模式。

Hermiticity 判定及紀錄：

- `PASS_STRUCTURE`：原生 ancilla-only PREP 與 unprepare 精確互逆，全域相位匹配，
  且原生 SELECT 的每個 ancilla basis branch 經精確 symbolic execution 證明
  無 leakage、為 Hermitian operator，因此 P†SP 為 Hermitian。
- `PASS_EXACT_COEFFICIENTS`：實際 block 的所有 Pauli 係數虛部都是精確零。
- `PASS_WITHIN_DELTA`：無法通過精確模式，但所有 Pauli 係數的虛部上界嚴格小於 delta。
  此結果標記為近似接受；JSON 仍保留原本的精確判定與非零係數證據。
- `FAIL_DELTA`：未通過精確模式，且最大虛部上界不小於 delta。
- `FAIL_NON_HERMITIAN` 或 `UNPROVEN`：delta=0 時沿用的精確模式結果。
  前者有虛部區間排除零；後者表示區間含零但沒有精確證據。

原生 SELECT 的多控制分解可能包含 decimal pi 旋轉。即使 proportionality 在 epsilon_BE
內通過，也可能不滿足精確 Hermiticity。新的容忍模式接受虛部上界低於 delta 的案例，
並以 `exact_hermiticity_status`、`hermiticity_status` 分別記錄精確結論與實際採用的判定。
這個 delta 門檻衡量 Pauli 係數，精確 Hermiticity 的 TASK 2 假設仍需另外處理。
生成時的小矩陣比較只檢查匯出與來源數值一致，不充當證明。
certificate 的範圍是 QASM 中精確寫下的小數與記錄的全域相位。

## 產物

- `benchmarks/pennylane/manifest.csv`：輸入索引與 metadata SHA-256。
- `benchmarks/pennylane/ising/ising_periodic_n*/`：原生 QASM、structure、來源／版本／hash metadata。
- `evaluation/pennylane/results/block_encoding_results.csv`：時間、記憶體、alpha、rigorous residual、
  delta、最大虛部上界、精確／近似 Hermiticity 與整體 verdict。
- `evaluation/pennylane/results/cases/*.json`：每列具體 block 係數、Hermiticity 證據與接受規則。
- `evaluation/pennylane/results/run_metadata.json`：命令、環境、manifest 與 verifier source hashes。

新執行的 runner 會把所有計入 source hashes 的 Python 檔案保存到輸出資料夾的
`source_snapshot/`，並在 metadata 的 `source_snapshots` 記錄對應位置。
`input_manifest.csv` 保存該次 manifest 的原始 bytes；其中的相對輸入路徑仍以
metadata 記錄的原始 `manifest_path` 所在目錄為基準。
README 的重跑範例指定 `/private/tmp` 輸出，以保留已保存的量測批次。

此變更的重跑結果另存於 `results_delta/`，保留先前的 `results/`。
執行中的 Python 程序要重新啟動，才會一致使用新的容忍模式和 CSV 欄位。

重現 n=2,3,4 的 delta 校準套件：

```bash
.venv/bin/python evaluation/pennylane/run_block_encoding.py \
  --manifest benchmarks/pennylane/delta_calibration_manifest.csv \
  --hermiticity-delta 1e-12 --block-encoding-epsilon 1e-8 \
  --precision 100 --output-dir /private/tmp/pennylane-ising-delta
```

runtime 包含 input validation、symbolic execution、alpha／residual 與 Hermiticity；
不包含生成、worker process startup 或最後 report 寫入。
`symbolic_runtime_sec`、`hermiticity_runtime_sec` 另列。每列 timeout 則記錄包含 process startup
的等待時間。ERROR／TIMEOUT 導致 runner exit code 1；測得語義 REJECT 的 exit code 是 0。

測試：

```bash
.venv/bin/python -m pytest tests/test_pennylane_block_encoding.py -q
PYTHONPATH=. ../cobble/.venv/bin/python -m pytest \
  tests/test_pennylane_block_encoding.py -q -k generation
```
