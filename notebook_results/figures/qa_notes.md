# 图件核验

- 绘图后端：Python / Matplotlib；宽度 182.9 mm。SVG 保留文本，PDF 嵌入 TrueType 字体，四份 PDF 的最小文本字号均为 8 pt。
- 数据来自当前四节点合成算例及 Gurobi 实际计算结果；未删除节点、走廊或有效前沿点。结果为确定性优化计算，无抽样误差条。
- 网络图直接读取 `corridors` 的端点、长度和现状标记，以及 `config.loads_kw` 和配变容量。包含全部 4 个节点、3 条现有线路、2 条候选走廊；布局坐标仅示意拓扑。
- 前沿图对应 `exact_frontier.csv`，切割算法恢复的前沿在内存中交叉核验；收敛图对应 `cold_start_cut_trace.csv`；快照图对应 `cut_snapshot_values.csv`。快照预算采样仅用于绘图。
- 已逐图检查节点、长度、图例、轴标签和曲线，无裁切或文字重叠。快照图用虚线叠在深色实线上表示两条结果重合。
- 源码预检 18 项通过、0 项失败。两项提示为未导出 TIFF、PNG 为 300 dpi：本次用途为 Notebook 显示及矢量导出，PNG 仅作预览，使用 SVG/PDF 获取矢量图。

| 图件 | 科学问题 | 视觉核验 |
|---|---|---|
| 00_network_topology | 当前节点、负荷与可用走廊如何连接 | 通过 |
| 01_planning_domain | 整数前沿与二维凸包有何差别 | 通过 |
| 02_cut_convergence | 固定预算的上下界是否收敛 | 通过 |
| 03_cut_snapshots | 累计电气割如何收紧整数主问题 | 通过 |
