"""MATPOWER case33bw 原始径向网架；三个实际节点的负荷独立变化。"""
from pathlib import Path  # 根据当前网架文件定位原始算例数据。
from copy import copy  # 实例化单个方案时复制母网，避免修改原始工况。
from functools import cached_property  # 不随查询变化的型号表和对照方案表只生成一次。
from itertools import product  # 仅用于小规模枚举核对，不进入紧凑求解器。
import re  # 从 MATPOWER 文本中提取数据表。

import numpy as np  # 保存节点、支路和型号参数数组。

from . import LineOptions, Project, RadialNetwork  # 复用统一的改造项目、型号表及径向网架结构。


class Case33(RadialNetwork):  # 在原始 33 节点径向网中加入逐线路改造选项。
    # 合成投资试验：费用为相对单位，不是 MATPOWER 提供的工程造价。
    projects = (Project("A", (2, 3), 2.), Project("B", (12, 13), 1.),  # 前两条候选支路及相对投资费用。
                Project("C", (23, 24), 1.), Project("D", (27, 28), 1.))  # 其余两条候选支路及相对投资费用。

    def __init__(self, load_nodes=(18, 25, 33)):  # 默认以实际节点 18、25、33 的负荷为三个独立坐标。
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
        design.x = np.asarray(choice, dtype=int)  # 四条线路分别取 0/1 型，与历史建设向量一致。
        design.cost = sum(o.cost[k] for o, k in zip(self.line_options, choice))  # 按选中的逐线路型号累加实际增量投资。
        return design  # 返回该次选型对应的固定网架。

    @cached_property  # 完整对照方案只生成一次，正式求解器不访问它。
    def designs(self):  # 为当前四条候选线路生成 16 方案独立基准。
        """仅供小算例枚举对照；紧凑 MILP/MISOCP 与联合割不读取此属性。"""
        designs = [self.design(choice) for choice in product(  # 逐线路型号的笛卡尔积仅用于这个小规模基准。
            *(range(len(o.cost)) for o in self.line_options))]  # 各线路的型号数直接读取唯一选项表。
        return tuple(sorted(designs, key=lambda d: (d.cost, tuple(d.x))))  # 按费用升序扫描可获得独立 AC 的最小可行投资。


network = Case33()  # 导出 Notebook 可直接加载的算例配置。
