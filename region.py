"""非负纯负荷径向网的单调区域分类；不跨建设方案取凸包。"""
import numpy as np  # 网格索引和成块分类使用同一坐标顺序。


def classify_orthant(states, index, status):  # 纯负荷单调域中，用一个已证点分类一块网格中心。
    """调用者须验证非负负荷、正阻抗和 vmax≥根电压；不跨方案构造凸包。"""
    if status == 1:  # 同一建设方案在更低负荷下仍可行。
        states[tuple(slice(0,int(i)+1) for i in index)] = 1  # 只覆盖坐标逐项不大于认证点的网格中心。
    elif status == -1:  # 若所有预算内方案在该点不可行，更高负荷也不可能可行。
        states[tuple(slice(int(i),None) for i in index)] = -1


def split_grid_box(lower, upper):  # 把含网格中心的闭索引盒沿最长轴分成两块。
    axis = int(np.argmax(upper-lower))  # 索引尺度对应相同的每轴网格分辨率。
    middle = (lower[axis]+upper[axis])//2  # 整数二分不会丢失或重复中心点。
    left, right = upper.copy(), lower.copy()  # 保留其他两个坐标的完整索引范围。
    left[axis], right[axis] = middle, middle+1  # 左右子盒不重叠，且并集等于父盒。
    return ((lower,left),(right,upper))
