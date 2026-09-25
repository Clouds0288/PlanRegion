# 项目清理记录

2026-09-25：在恢复对偶主线后，将五节点勘察中仍需使用的函数整理为正式模块，清除临时实验入口和历史试跑目录。

- `survey.py`：按需求域、共享割、道路评分、逐轮勘察和结果导出。参数直接来自正式网架，不读取旧实验结果或测试夹具。
- `plot.py`：勘察绘图统一由 `render_survey(result, output)` 完成，原 `survey_plot.py` 已合并删除；保留三维结果和离线回放。
- `Network/`：保留五节点、四节点、Case33 和江口网架及其源数据。
- `tests/`：保留对偶、物理、几何、入口和界面回归。两个仍有测试依赖的比较函数提取到 `planning_checks.py`，其余一次性 benchmark、probe、audit、repair 和独立实验脚本归档。
- `docs/notation-history.md`：完整保存清理前符号契约；活动定义及接口退役见 `notation.md` 第 18 节。未改数学量名称、单位、切片或评分公式。
- `results/`：保留已审核的 Case33 latest 数据及五节点对偶复现结果。失败试跑、旧枚举与其他探索性实验移出工作区。

清理前资料保存在项目同级目录：

`../PlanRegion-cleanup-archive/20260925-174418/`

其中 `before_cleanup.zip` 为可恢复归档，`manifest.json` 记录原路径、大小和 SHA-256，`working_changes.patch` 保存未提交改动。归档包含 1396 个文件，约 1098.95 MB 原始内容，压缩后约 254.64 MB；已校验 ZIP 完整性。保留中的 Case33 latest 数据、Git 元数据、运行环境和工具缓存不重复放入归档。

GitNexus 的依赖与影响分析在清理前执行；没有将 UNKNOWN 当作无依赖。主入口、Case33 两阶段实验及共享校验文件的 HIGH 风险、面积增量函数的 CRITICAL 风险均已说明。清理顺序为提取必要实现、更新调用方、验证活动依赖，再移除旧目录。

实际采用可恢复的移出归档方式：96 个文件或目录目标、1302 个文件（其中 158 个 Python 文件），约 1096.5 MB，已移出项目。仅供旧 benchmark 使用的 `plot.save_benchmark_report` 一并退役，旧实现包含在清理前源码归档中。经用户确认，结束了 9 月 23 日启动且占用历史四节点结果的旧进程。

验证：完整默认回归 114 项通过、234 项子测试通过、2 项按配置跳过；包含五节点脱离旧实验目录后的完整重算。迁移后两幅五节点 PNG 与保留的已审核图逐像素一致。数学符号检查另行通过。

重新复现：`python survey.py --threads 1`。完整回归：`python -m pytest tests -q`。
