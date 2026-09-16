"""从空内核完整执行 notebook，并保存全部输出。"""
from pathlib import Path
import asyncio

import nbformat
from nbclient import NotebookClient


path = Path(__file__).with_name("planning_domain_demo.ipynb")
notebook = nbformat.read(path, as_version=4)
client = NotebookClient(
    notebook,
    timeout=600,
    kernel_name=notebook.metadata.kernelspec.name,
    resources={"metadata": {"path": str(path.parent)}},
    extra_arguments=["--HistoryManager.hist_file=:memory:"],
)
with asyncio.Runner(loop_factory=asyncio.SelectorEventLoop) as runner:
    runner.run(client.async_execute())
nbformat.write(notebook, path)
print(f"已完成 {sum(c.cell_type == 'code' for c in notebook.cells)} 个代码单元：{path.name}")
