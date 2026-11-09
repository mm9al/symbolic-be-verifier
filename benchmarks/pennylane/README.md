# Native PennyLane block encodings

`manifest.csv` 索引原生 PennyLane `PrepSelPrep` periodic-Ising 電路。
每組只保存 QASM、PREP／SELECT／unprepare 結構與來源 metadata；實驗結果在
`evaluation/pennylane/results`。

生成、TASK 1 接受條件與重跑指令見
[evaluation/pennylane/README.md](../../evaluation/pennylane/README.md)。
主 manifest 的尺寸由最近一次 generator 的 `--sizes` 控制。
`delta_calibration_manifest.csv` 固定列出 n=2,3,4，供精確／容忍 Hermiticity 模式比較。
`all_ising_manifest.csv` 索引目前已保存的完整集合：
n=2..10,16,17,32,33,64,65,128。重跑已保存的實驗時明確指定這兩份固定索引，
避免使用會隨最近一次生成改變的主 manifest。
