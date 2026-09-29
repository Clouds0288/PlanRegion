"""Case33 固定功率因数正负接入：二维60秒、三维300秒；独立扫描和完整回放。

python experiments/case33_signed.py
python experiments/case33_signed.py --dimension 3
python experiments/case33_signed.py --replay results/case33_signed/mode_1/18_25_30/recording.json.gz
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Network.case33bw import Case33
from experiments.fourbus_signed import main

# 1. 全部符号区共用构域时限；未选中的节点保留原始背景负荷。
LOAD_NODES = (18, 25, 30)
TIME_LIMITS = {2: 60., 3: 300.}
BUDGET = 7.
OUTPUT = ROOT/'results'/'case33_signed'

# 2. 与 FourBus 共用构域、扫描、缓存和逐帧回放，不复制算法。
if __name__ == '__main__':
    main(network_type=Case33, load_nodes=LOAD_NODES, budget=BUDGET,
         time_limits=TIME_LIMITS, output=OUTPUT)
