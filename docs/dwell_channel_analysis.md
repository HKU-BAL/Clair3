# 全对齐张量信号通道差异分析

## 背景

在为 Clair3 的全对齐张量添加来自 BAM `mv` 标签的信号长度（dwell time）通道时，我们同时维护了 Python 实现（`CreateTensorFullAlignment.py`）与新的 C/CFFI 实现（`clair3_full_alignment_dwell.c` + `CreateTensorFullAlignmentFromCffiWithDwell.py`）。为了让两端生成的张量完全一致，需要确保 **mv 解析、插入信号累积以及张量写入** 的逻辑完全对齐。

## 差异原因

初期对比发现所有差异都集中在新的第 9 个信号通道（channel 8，0-based）上，主要来源如下：

1. **归一化与截断**：C 端最初对信号长度做了缩放（除以 10）与 `int8` 截断，而 Python 端保留原始计数，导致数值系数不一致。
2. **插入片段处理**：Python 在遍历 CIGAR 时，会把一整个插入片段的 dwell time 累积到插入位置对应的参考坐标；C 端最初把插入信号存放在独立字段，张量填充时又覆盖掉匹配/错配位置，逻辑不一致。
3. **mv 标签解析**：Python 端接收到的 `mv_tag` 已经去掉了 stride 信息（`mv:B:c,<0/1序列>`），而 C 端仍试图读取并展开 stride，导致每个碱基的 dwell 计数错位。
4. **读排序差异**：在修复前三项后，仍有少量差异。这是因为 Python 端在 mpileup 结果基础上排序读，而 C 端在遍历读覆盖面时采用了不同的排序方式。对张量按 depth 维度逐行排序后，两端数据完全一致。

## 关键修改

1. **信号长度字段**
   - `Pos_info.signal_length` 类型改为 `int32_t`，并去掉 `ins_signal_length` 字段。
   - 插入事件的 dwell time 直接累加到 `signal_length` 中，保持与 Python 相同的“插入挂在参考坐标”的策略。

2. **mv 标签解析**
   - `compute_signal_lengths_from_mv_tag` 仅使用 `mv` payload，不再做 stride 展开。
   - 遇到 `1` 写当前 dwell 计数，遇到 `0` 继续累加，并在负链时整体翻转数组，与 Python 的 `compute_dwell_times` 逐元素一致。

3. **张量写入**
   - 取消 `encode_signal_value` 缩放/截断逻辑，直接把 `signal_length` 按 `int8_t` 范围截断（0~127）写入第 9 通道。
   - 插入事件写在同一位置，不额外产生独立通道。

4. **CFFI 包装**
   - 新增 `CreateTensorFullAlignmentFromCffiWithDwell.py`，根据 `--enable_dwell_time` 参数自动设置通道数并调用新的 C 入口。
   - `build_dwell.py` 链接新的 `clair3_full_alignment_dwell.c`，生成 `libclair3_dwell.so` 并拷贝到仓库根目录。

## 验证流程

1. 重新运行 Python 与 C 版本的全对齐张量生成脚本（均开启 `--enable_dwell_time`），并使用 `np.load`/文本解析方式读取结果。
2. 比较 `abs(c_tensor - py_tensor)`，确认只有 dwell 通道存在差异，再对张量在 depth 维度进行排序，确保排序后的结果完全一致，以排除读顺序差异。
3. 针对个别读，利用 `pysam` 调用 `compute_dwell_times` 与 C 端 `compute_signal_lengths_from_mv_tag` 的副本进行逐读对比，确认两者输出一致。

## 剩余差异与结论

最终仅剩的 3 处差异来自读排序。对张量行进行排序（或按 `(haplotype, read_name)` 等键规整顺序）后，Python 与 C 的张量完全一致，说明逻辑层面的差异已被消除。

## 建议

1. **测试脚本**：对比张量时建议加入排序或根据读名对齐的逻辑，以避免不同实现之间的读顺序差异造成假阳性。
2. **代码维护**：在文档/注释中记录 `mv` 标签无需 stride 展开的事实，避免后续开发者引用旧逻辑。
3. **扩展**：如需在更多流程启用 dwell 通道，可复用当前 CFFI 封装；若要进一步加速，可考虑在 C 端直接按读名排序。

本次修复后，信号通道已可稳定用于训练与推理，且 Python/C 实现可相互校验。

