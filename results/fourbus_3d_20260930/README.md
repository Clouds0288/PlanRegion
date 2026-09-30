# FourBus 3D 正负功率结果归档

本目录是 2026-09-30 正负功率（`mode=1`）FourBus 三维实验的正式索引。两种方法均使用节点 `(1,2,3)`、预算 20000、20 个求解线程、总构域时限 300 秒和 `tau=0.005`。

## 共同扫描参考

- 文件：`../scans/fourbus_signed/mode_1/1_2_3/pf_0.95_1_+1_budget_20000.npz`
- SHA-256：`d2f247758b001d9e9c6e8354ad9a9f0defc85d6dce6348bb1207f3b6c820160f`
- 网格：`80×80×80`，共 512000 个格点，其中 184784 个可行点
- 范围：`[-170,-135,-135]` 至 `[158,128,127]` kW

扫描只用于两种结果的同网格评价，不参与构域或终止判定。

## 当前主线

- 求解结果：`../support_face_test/fourbus_3d_signed_trial/mainline_1/result.json`
- 原生回放：`../support_face_test/fourbus_3d_signed_trial/mainline_1/monitor.json.gz`
- 配置、源码哈希和对比指标：`../support_face_test/fourbus_3d_signed_trial/protocol.json`、`comparison.json`、`report.md`
- 构域耗时：271.02 秒（对比程序记录的墙钟时间为 272.62 秒）
- 覆盖认证：1/8 个功率符号分区完成，整体因时限未认证
- 扫描指标：内域遗漏率 5.66283%，误收率 0%；外域遗漏率 0%，误收率 44.48830%

项目主线另有直接运行记录：`../mainline/mode_1/fourbus_1_2_3.json.gz`。

## 新支持面测试方案（完整 physical）

- 求解结果：`../support_face_test/fourbus_3d_physical_verified/support_1/result.json`
- 逐分区证书：同目录的 `partition_*.json`
- 原生回放：同目录的 `monitor.json.gz`
- 固定源码、配置、见证审计和指标：`../support_face_test/fourbus_3d_physical_verified/`
- 构域耗时：209.49 秒（墙钟时间 210.67 秒）
- 覆盖认证：5/8 个功率符号分区完成，整体因三个分区的全局 physical 查询超时而未认证
- 扫描指标：内域遗漏率 0.07468%，误收率 0%；外域遗漏率 0%，误收率 30.39154%

新方案的 1439 个支持点和全局见证全部通过原物理约束复核，最大违反量为 `1.2616e-9`。严格结论与适用边界以两份结果目录中的报告和 JSON 为准。
