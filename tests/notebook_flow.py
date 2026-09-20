"""测试直接执行主 Notebook 的函数单元，避免另写一套实验算法。"""
from functools import lru_cache  # 一个测试进程只加载一次 Notebook 函数。
import nbformat  # 直接读取主入口，避免维护第二套流程。


@lru_cache  # 缓存已执行的定义单元命名空间。
def workflow():  # 为测试加载真实 Notebook 中的求解及验证函数。
    notebook = nbformat.read('main.ipynb', as_version=4)  # 按标准 Notebook 格式读取源单元。
    namespace = {}  # 所有定义单元共享同一独立命名空间。
    tags = {'imports', 'linear-flow', 'cut-flow', 'method-flow', 'evaluation-flow', 'joint-flow'}  # 仅加载依赖与函数，不启动完整实验或绘图。
    for cell in notebook.cells:  # 保持 Notebook 中定义单元的原有顺序。
        if tags.intersection(cell.metadata.get('tags', [])):  # 依据显式标签选择可测试的流程单元。
            exec(compile(cell.source, 'main.ipynb', 'exec'), namespace)  # 测试直接执行正式流程代码，报错标明主入口来源。
    return namespace  # 返回函数和依赖供回归测试调用。
