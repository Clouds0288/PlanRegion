"""测试直接执行主 Notebook 的函数单元，避免另写一套实验算法。"""
from functools import lru_cache
import nbformat


@lru_cache
def workflow():
    notebook = nbformat.read('main.ipynb', as_version=4)
    namespace = {}
    tags = {'imports', 'linear-flow', 'cut-flow', 'method-flow', 'evaluation-flow', 'joint-flow'}
    for cell in notebook.cells:
        if tags.intersection(cell.metadata.get('tags', [])):
            exec(compile(cell.source, 'main.ipynb', 'exec'), namespace)
    return namespace
