# 数学符号与变量契约

- 修改数学代码前必须阅读 [docs/notation.md](docs/notation.md)。它约定数学符号与代码名、单位、形状、索引，以及持久化与结果格式。
- 登记只覆盖三类名字：代表数学量或数值设置的名字（类字段与属性、容差与默认常量）；持久化与结果格式的键（回放记录、扫描缓存、导出文件，以及交给其他模块的结果）；公开的函数与类。函数内的局部量、形参与下划线开头的私有名是实现细节，不登记，`tests.test_notation` 会拒绝这类登记。
- 已登记的名字在语义未变时保留；改名、改单位、符号方向、索引、切片或格式须在 notation.md 写显式迁移，同步调用方、文档与测试；不能删登记项来绕过检查。新增数学量、持久化键或公开接口先登记。
- 实现细节可以自由重构：拆函数、改局部量或形参名、调整私有辅助函数，不需要登记，也不需要写迁移。
- 改名作为独立事项说明理由，遵守下方 GitNexus impact / rename 规则。rename 工具不可用时，先用 impact 与 context 列出全部调用方和引用，逐处手工修改，再用文本搜索确认没有遗漏，最后运行 detect-changes 与测试；不做盲目的全局查找替换。历史源码快照和结果不做机械改名。
- 交付前运行 `python -m unittest tests.test_notation -v`；算法变更另运行相关数值回归。名称检查不替代数学正确性验证。

# 代码风格：精炼整洁，不过度防御

- 写科研代码追求精炼整洁：代码直接表达数学与算法，主线一眼可读。
- 只为真实会发生、且有明确处理方式的情况写分支；不为不可能出现的输入加缺省值、类型或形状判断、try/except 兜底。
- 失败要显式：求解或几何失败直接抛出，或抛专门的异常由调用方按“未决”处理；不静默吞掉异常，不用缺省值掩盖。
- 外部输入（命令行参数、读入的文件）在边界处校验一次；内部调用之间不重复校验。
- 不保留没有调用方的参数、模式、返回字段和兼容分支；涉及已登记名字时按上面的契约写迁移。

<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **PlanRegion** (216 symbols, 261 relationships, 5 execution flows).

> Index stale? Run `node .gitnexus/run.cjs analyze --index-only` from the project root — it auto-selects an available runner. No `.gitnexus/run.cjs` yet? Bootstrap with `npx`, `bunx`, or `pnpm dlx` — e.g. `bunx gitnexus@latest analyze` (npm 11 npx crash; #1939).

## Always Do

- **MUST run impact before editing.** Use `impact({target: "symbolName", direction: "upstream"})` or `node .gitnexus/run.cjs impact "symbolName" --direction upstream --repo .`; report callers, processes, and risk. Never substitute grep for graph analysis.
- **MUST analyze graph changes before committing.** Use `detect_changes({scope: "all"})` (MCP) or `node .gitnexus/run.cjs detect-changes --scope all --repo .` (CLI fallback). `partial: true` or `truncated: true` is not a clean check — a zero means unseen, not unaffected; re-run it. For regression review: `detect_changes({scope: "compare", base_ref: "main"})` or `node .gitnexus/run.cjs detect-changes --scope compare --base-ref "main" --repo .`.
- MUST warn on HIGH/CRITICAL `risk` pre-edit; never use `riskSharedAxes` to waive a HIGH/CRITICAL `risk` warning. Compare File/symbol: MCP File omits axes; Graph-RAG expands File.
- **MUST treat `risk: UNKNOWN` as unresolved, not as low.** An empty caller set is not evidence the symbol is unused — it can also mean the callers are not resolvable by the index (plain-object property access, dynamic dispatch, cross-language calls). `impact` pairs `UNKNOWN` with a `riskNote` saying so. Confirm with a text search before treating the symbol as safe to change or delete; do not proceed on the strength of a zero.
- **MUST use `query({search_query: "concept"})` for concepts/flows, `context({name: "symbolName"})` for a named symbol, or `impact` for blast radius, on read-only callers, dependencies, imports, or execution flow.** Graph first; text search only for empty/`UNKNOWN`/literals.
- For security review, `explain({target: "fileOrSymbol"})` lists taint findings (source→sink flows; needs `analyze --pdg`).

## Never Do

- NEVER edit a function, class, or method before MCP/CLI impact analysis.
- NEVER ignore HIGH or CRITICAL risk warnings from impact analysis, and never read `UNKNOWN` as an all-clear — it means the walk could not answer, which is the one verdict that requires confirming by other means.
- NEVER rename symbols with find-and-replace — use `rename` which understands the call graph.
- NEVER commit before MCP/CLI graph change analysis.

## Resources

| Resource | Use for |
| --- | --- |
| `gitnexus://repo/PlanRegion/context` | Codebase overview, check index freshness |
| `gitnexus://repo/PlanRegion/clusters` | All functional areas |
| `gitnexus://repo/PlanRegion/processes` | All execution flows |
| `gitnexus://repo/PlanRegion/process/{name}` | Step-by-step execution trace |

## CLI

| Task | Read this skill file |
| --- | --- |
| Understand architecture / "How does X work?" | `.claude/skills/gitnexus-exploring/SKILL.md` |
| Blast radius / "What breaks if I change X?" | `.claude/skills/gitnexus-impact-analysis/SKILL.md` |
| Trace bugs / "Why is X failing?" | `.claude/skills/gitnexus-debugging/SKILL.md` |
| Rename / extract / split / refactor | `.claude/skills/gitnexus-refactoring/SKILL.md` |
| Tools, resources, schema reference | `.claude/skills/gitnexus-guide/SKILL.md` |
| Index, status, clean, wiki CLI commands | `.claude/skills/gitnexus-cli/SKILL.md` |

<!-- gitnexus:end -->
