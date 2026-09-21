"""兼容既有测试接口；算法的唯一来源是 main.py。"""
from functools import lru_cache
import main


@lru_cache  # 缓存定义命名空间，不运行完整实验。
def workflow():  # 提供正式 MP/SP、AC 和构域函数供测试直接调用。
    return vars(main)  # 保留 patch.dict 对函数全局命名空间的测试替换。
