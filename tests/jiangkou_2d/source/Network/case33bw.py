"""MATPOWER case33bw 原始径向网架；三个实际节点的负荷独立变化。"""
from pathlib import Path  # 根据当前网架文件定位原始算例数据。
from copy import copy  # 实例化单个方案时复制母网，避免修改原始工况。
from functools import cached_property  # 不随查询变化的逐线路型号表只生成一次。
import re  # 从 MATPOWER 文本中提取数据表。

import numpy as np  # 保存节点、支路和型号参数数组。

from . import LineOptions, Project, RadialNetwork  # 复用统一的改造项目、型号表及径向网架结构。


class Case33(RadialNetwork):  # 在原始 33 节点径向网中加入逐线路改造选项。
    # 合成投资试验：费用为相对单位，不是 MATPOWER 提供的工程造价。
    projects = (  # 前 4、8、16、32 项形成嵌套候选集；费用与线路只在此定义。
        Project("A", (2, 3), 2.),     # 原候选：公共上游线路。
        Project("B", (12, 13), 1.),   # 原候选：节点 18 所在支路。
        Project("C", (23, 24), 1.),   # 原候选：节点 25 所在支路。
        Project("D", (27, 28), 1.),   # 原候选：节点 33 所在支路。
        Project("E", (5, 6), 1.),     # 新候选：下游分支的上游线路。
        Project("F", (13, 14), 1.),   # 新候选：进一步改善节点 18 方向。
        Project("G", (24, 25), 1.),   # 新候选：节点 25 的入线。
        Project("H", (28, 29), 1.),   # 新候选：进一步改善节点 33 方向。
        Project("I", (3, 4), 1.),     # 十六线路新增：节点 18、33 方向的公共上游。
        Project("J", (6, 7), 1.),     # 十六线路新增：节点 18 所在分支的入口。
        Project("K", (9, 10), 1.),    # 十六线路新增：节点 18 方向的中段。
        Project("L", (16, 17), 1.),   # 十六线路新增：节点 18 方向的末段。
        Project("M", (3, 23), 1.),    # 十六线路新增：节点 25 所在分支的入口。
        Project("N", (6, 26), 1.),    # 十六线路新增：节点 33 所在分支的入口。
        Project("O", (29, 30), 1.),   # 十六线路新增：节点 33 方向的中下游。
        Project("P", (31, 32), 1.),   # 十六线路新增：节点 33 方向的末段。
        Project("Q", (1, 2), 1.),
        Project("R", (4, 5), 1.),
        Project("S", (7, 8), 1.),
        Project("T", (8, 9), 1.),
        Project("U", (10, 11), 1.),
        Project("V", (11, 12), 1.),
        Project("W", (14, 15), 1.),
        Project("X", (15, 16), 1.),
        Project("Y", (17, 18), 1.),
        Project("Z", (2, 19), 1.),
        Project("AA", (19, 20), 1.),
        Project("AB", (20, 21), 1.),
        Project("AC", (21, 22), 1.),
        Project("AD", (26, 27), 1.),
        Project("AE", (30, 31), 1.),
        Project("AF", (32, 33), 1.),
    )

    def __init__(self, load_nodes=(18, 25, 33), candidate_count=4):  # 默认保留四线路基准，扩展时显式给出候选数量。
        if not 1 <= candidate_count <= len(type(self).projects):  # 数量由唯一项目表决定，后续扩展无需修改限制。
            raise ValueError(f'candidate_count must be between 1 and {len(type(self).projects)}')  # 禁止切片静默截短候选集。
        self.projects = type(self).projects[:candidate_count]  # 候选参数只记录一次，实例选择需要的前缀。
        path = Path(__file__).parent/"data"/"case33bw.m"  # 原始 MATPOWER 文件是基础工况的唯一来源。
        source = path.read_text(encoding="utf-8")  # 读取原始文本，不另存一份参数副本。
        tables = {}  # 保存从文本解析出的三张标准数据表。
        for name in ("bus", "gen", "branch"):  # 依次读取节点、电源和支路表。
            block = re.search(rf"mpc\.{name}\s*=\s*\[(.*?)\];", source, re.S)[1]  # 提取对应 mpc 字段的方括号内容。
            block = re.sub(r"%[^\n]*", "", block)  # 去除 MATLAB 行注释后再解析数值。
            tables[name] = np.array([[float(x) for x in row.split()]  # 按空白拆分每行字段，并转为浮点数。
                                     for row in block.split(";") if row.strip()])  # 按分号拆分标准数据行，跳过末尾空段。
        bus, gen, branch = (tables[k] for k in ("bus", "gen", "branch"))  # 按固定字段名称取得三张表。
        base_mva = float(re.search(r"mpc.baseMVA\s*=\s*([\d.]+)", source)[1])  # 读取原始功率基准，不能把它当作变压器容量。
        root = int(bus[bus[:, 1] == 3, 0][0])  # MATPOWER 类型 3 的节点是电源根节点。
        rows = bus[bus[:, 0] != root]  # 其余节点全部进入非根节点表。
        nodes = tuple(map(int, rows[:, 0]))  # 保留真实节点编号，供显示和独立交叉核验。
        selected = [nodes.index(i) for i in load_nodes]  # 独立变化的真实负荷节点对应的数组位置。
        active = branch[branch[:, 10] == 1]  # 五条常开联络线保持断开。
        zbase = bus[0, 9]**2/base_mva  # 原始线路阻抗单位是 Ω，用 Ubase²/Sbase 换成标幺值。
        power_limit = gen[0, 8]*1000-rows[:, 2].sum()+rows[selected, 2].sum()  # 电源有功上限减固定背景负荷，得到独立节点无损总量外界。
        # rateA=0：原始数据未启用线路热限；baseMVA 不是配变容量。
        super().__init__(  # 将原始数据转换到所有模型共享的径向接口。
            "case33bw", root, nodes, active[:, :2].astype(int), active[:, 2]/zbase,  # 在运支路电阻由 Ω 换算为标幺值。
            active[:, 3]/zbase, base_mva*1000, bus[0, 9], rows[:, 2], rows[:, 3],  # 同样转换电抗；原始节点负荷继续保留 kW/kvar。
            tuple(load_nodes), rows[selected, 3]/rows[selected, 2], rows[:, 12]**2,  # 由原工况确定独立节点 Q/P 比，并将电压下限平方。
            rows[:, 11]**2, power_limit, source_pmax=gen[0, 8]/base_mva,  # 电压幅值限值先平方；Pmax 从 MW 转标幺。
            source_qmax=gen[0, 3]/base_mva,  # Qmax 从 Mvar 转标幺，保留原始电源运行限制。
            sources=("Network/__init__.py", "Network/case33bw.py", "Network/data/case33bw.m"))  # 实验记录使用这些源文件的指纹追踪数据版本。

    @cached_property  # 逐线路型号表由母网配置派生一次。
    def line_options(self):  # 向紧凑模型提供逐支路选型输入。
        """紧凑模型只读取逐线路选项，不生成完整建设组合。"""
        options = []  # 收集候选支路各自的型号表。
        for project in self.projects:  # 投资项目的支路与费用只在 projects 中定义。
            e = self.nodes.index(project.branch[1])  # 入边以受端节点索引标识，拓扑保持不变。
            options.append(LineOptions(e, (self.r[e], self.r[e]/2),  # 两种电阻分别对应保持原状与并联一回。
                                       (self.reactance[e], self.reactance[e]/2),  # 同样按两回相同线路并联的规则计算电抗。
                                       (0., project.cost)))  # 同走廊并联一回，R/X 减半；费用仍由 projects 唯一定义。
        return tuple(options)  # 返回可复用的逐线路型号配置。

    def design(self, choice):  # 仅按给定选型创建一个物理网架实例。
        """按一次求解给出的线路型号实例化网架；不遍历其他组合。"""
        design = copy(self)  # 固定拓扑与背景数据可以共享，只复制实例外壳。
        design.r, design.reactance = self.r.copy(), self.reactance.copy()  # 阻抗单独复制，因为当前建设方案会改变这些值。
        for k, options in zip(choice, self.line_options):  # 按每条候选线路选中的型号更新物理参数。
            design.r[options.branch] = options.r[k]  # 替换该支路的标幺电阻。
            design.reactance[options.branch] = options.reactance[k]  # 替换该支路的标幺电抗。
        design.x = np.asarray(choice, dtype=int)  # 每条候选线路的型号编号，长度随候选数量变化。
        design.cost = sum(o.cost[k] for o, k in zip(self.line_options, choice))  # 按选中的逐线路型号累加实际增量投资。
        return design  # 返回该次选型对应的固定网架。


network = Case33()  # 导出 Notebook 可直接加载的算例配置。
