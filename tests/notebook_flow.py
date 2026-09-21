"""测试读取正式 Notebook 函数；审核和看门狗不进入主线。"""
from functools import lru_cache  # 同一测试进程只加载一次函数。
import nbformat  # 以 Notebook 为流程的唯一来源。


@lru_cache  # 缓存定义命名空间，不运行完整实验。
def workflow():  # 提供正式 MP/SP、AC 和构域函数供测试直接调用。
    namespace = {}  # 所有函数共享导入环境。
    for cell in nbformat.read('main.ipynb',as_version=4).cells:  # 按定义顺序读取单元。
        tags = cell.metadata.get('tags',[])  # Notebook 只标记导入和三个流程单元。
        if cell.cell_type=='code' and ('imports' in tags or any(t.endswith('-flow') for t in tags)):  # 跳过配置、执行及绘图。
            exec(compile(cell.source,'main.ipynb','exec'),namespace)  # 测试直接执行正式代码。
    return namespace  # 不复制一套算法。
