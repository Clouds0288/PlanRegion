"""从唯一 history 重建内外域；变化后的显示网格缓存一次，播放不重新求解。"""
import json
from pathlib import Path

import numpy as np
from IPython.display import IFrame, display
from plotly.offline import get_plotlyjs

from region import add_certificate, simplex, update_outer
from plot import PLOT_SCRIPT, certified_frame, frame_data


class Replay:
    """同一网架的多预算回放；每档保留独立历史，网架参数只打包一次。"""

    def __init__(self, processes):
        self.processes = processes
        self.states = []
        for process in processes:
            polytopes = [simplex(process.network) for _ in process.designs]
            self.states.append(dict(polytopes=polytopes, certified={},
                                    frames=[frame_data(process.costs, process.designs, polytopes)],
                                    inner_frames=[certified_frame({})], timeline=[[0, 0]]))

    def data(self):
        scenarios = []
        for process, state in zip(self.processes, self.states):
            for row in process.history[len(state["timeline"])-1:]:
                updated = update_outer(process.costs, process.designs, state["polytopes"],
                                       row, process.queries[row["query"]])
                if any(not np.array_equal(old, new) for old, new in zip(state["polytopes"], updated)):
                    state["frames"].append(frame_data(process.costs, process.designs, updated, row["cut"]))
                state["polytopes"] = updated
                if add_certificate(state["certified"], row):
                    frame = certified_frame(state["certified"])
                    if frame != state["inner_frames"][-1]:
                        state["inner_frames"].append(frame)
                state["timeline"].append([len(state["frames"])-1, len(state["inner_frames"])-1])
            scenarios.append(dict(budget=process.budget if np.isfinite(process.budget) else None,
                                  complete=process.finished, queries=process.queries, history=process.history,
                                  frames=state["frames"], inner_frames=state["inner_frames"],
                                  timeline=state["timeline"]))
        network = self.processes[0].network
        return dict(network=network.name, limit=network.power_limit, cost=network.cost,
                    load_nodes=network.load_nodes, scenarios=scenarios)

    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(self.data(), ensure_ascii=False, separators=(",", ":"),
                             allow_nan=False, default=lambda value: value.tolist())
        html = (PAGE.replace("__LIBRARY__", get_plotlyjs())
                .replace("__PLOT__", PLOT_SCRIPT).replace("__DATA__", payload.replace("</", "<\\/")))
        path.write_text(html, encoding="utf-8")
        return path

    def show(self, path):
        path = self.save(path)
        # Notebook 只引用这一份 HTML，不再把相同历史嵌入输出。
        display(IFrame(path.as_posix(), width="100%", height=1050))


PAGE = r"""<!doctype html>
<html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>节点自由分配规划 · 预算与切割回放</title>
<style>
body{font:14px system-ui,'Microsoft YaHei',sans-serif;max-width:1100px;margin:20px auto;padding:0 18px;color:#253b49;background:#fff}
h1{font-size:22px;font-weight:600;margin-bottom:8px}.controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:12px 0}
button,select,input{font:inherit}button,select{border:1px solid #bac6cd;border-radius:4px;padding:6px 10px;background:white;color:inherit}
button{cursor:pointer}button:disabled{opacity:.4}input[type=range]{flex:1;min-width:150px;accent-color:#48799a}
#round{width:65px}#frontier{height:290px}#region{height:480px}.caption{color:#526571;margin:4px 0 12px}
#status{min-height:22px}button:focus-visible,input:focus-visible,select:focus-visible{outline:2px solid #48799a}
</style><script>__LIBRARY__</script></head><body>
<h1>节点自由分配规划 · 预算与切割回放</h1><p id="parameters"></p>
<div class="controls"><label>预算 <select id="budget" aria-label="建设预算"></select></label>
<label>回放 <select id="mode" aria-label="步进方式"><option value="budgets">预算演变</option><option value="rounds">完整切割过程</option><option value="cuts">仅新增割</option></select></label></div>
<div class="controls"><button id="first">初始</button><button id="previous">上一轮</button>
<button id="play">播放</button><button id="next">下一轮</button><button id="last">最后</button>
<select id="speed" aria-label="播放间隔"><option value="100">0.1 秒</option><option value="250">0.25 秒</option><option value="1000">1 秒</option><option value="2000" selected>2 秒</option></select>
<input id="slider" type="range" min="0" step="1" aria-label="回放进度">
<label id="round-label">第 <input id="round" type="number" min="0" step="1" aria-label="指定轮次"> 轮</label></div>
<p id="status" role="status"></p><div id="frontier"></div>
<p class="caption">图1｜当前预算上限内的总负荷—最低预算前沿；橙点已通过 SP 认证。</p>
<div class="controls"><label><input id="outer" type="checkbox">显示候选外包络（可能含不可行点）</label>
<label><input id="before" type="checkbox">对比前一候选包络</label>
<span>拖动旋转，滚轮缩放；预算切换保持坐标尺度和视角。</span></div>
<div id="region"></div><p class="caption">图2｜三个轴为节点 1、2、3 的负荷，各节点可从零自由分配。蓝色为预算内各建设方案已认证可调度域的并集，保留非凸边界，不对所有方案取整体凸包。预算演变展示每档预算的最终进度；完整切割过程可从第0轮回放。未认证完成时蓝色仅为内近似。候选点判定针对本轮方案；结论针对当前线性电气模型。</p>
<script id="experiment" type="application/json">__DATA__</script><script>
const experiment=JSON.parse(document.getElementById('experiment').textContent);
__PLOT__
const el=id=>document.getElementById(id), scenarios=experiment.scenarios;
const budgetLabel=value=>value===null?'无限预算':value.toLocaleString()+' 元';
const budgetMode=()=>el('mode').value==='budgets';
const counts=scenarios.map(s=>{const out=[0];for(const row of s.history)out.push(out.at(-1)+(row.cut!==null?1:0));return out;});
const cutSteps=scenarios.map(s=>Array.from(new Set([0,...s.history.flatMap((r,i)=>r.cut!==null?[i+1]:[]),s.history.length])));
for(const [index,scenario] of scenarios.entries()){
 const option=document.createElement('option');option.value=index;option.textContent=budgetLabel(scenario.budget);el('budget').appendChild(option);
}
let selected=0,current=scenarios[0].history.length,timer=null;
const lastPosition=()=>budgetMode()?scenarios.length-1:scenarios[selected].history.length;
const position=()=>budgetMode()?selected:current;
function move(value){
 if(budgetMode()){selected=value;current=scenarios[selected].history.length;}else current=value;
 render();
}
function render(){
 const data=scenarios[selected], history=data.history;
 const row=history[current-1], count=counts[selected][current], points=[];
 for(const r of history.slice(0,current)){
   const query=data.queries[r.query];
   if(query.mode==='MP1-total' && r.eta!==null && r.cut===null && r.x!==null)
     points.push([query.target,experiment.cost.reduce((s,c,i)=>s+c*r.x[i],0)]);
 }
 const [outerIndex,innerIndex]=data.timeline[current];
 drawFrame(data.frames[outerIndex],data.frames[Math.max(0,outerIndex-1)],data.inner_frames[innerIndex],budgetMode()?null:row,points);
 el('parameters').textContent=`网架：${experiment.network} ｜ 建设预算：${budgetLabel(data.budget)} ｜ 节点负荷自由分配`;
 el('budget').value=selected;el('slider').max=lastPosition();el('slider').value=position();
 el('round').max=history.length;el('round').value=current;el('round-label').hidden=budgetMode();
 el('first').textContent=budgetMode()?'最低预算':'初始';el('previous').textContent=budgetMode()?'上一预算':'上一轮';
 el('next').textContent=budgetMode()?'下一预算':'下一轮';el('last').textContent=budgetMode()?'最高预算':'最后';
 const done=current===history.length && data.complete;
 const status=done?'完整规划—可调度域已认证':current===0?'初始候选域，尚无认证区域':row.x===null?'当前查询不可行':row.cut!==null?'新增可行性割，蓝色仍为已认证内域':'本轮方案通过 SP，完整区域仍待认证';
 el('status').textContent=`${budgetMode()?`预算 ${selected+1}/${scenarios.length} · `:''}${current} / ${history.length} 轮 · ${count} 条割 · ${status}`;
 el('previous').disabled=position()===0;el('next').disabled=position()===lastPosition();
}
function pause(){clearInterval(timer);timer=null;el('play').textContent='播放';}
function step(direction){
 const sequence=el('mode').value==='cuts'?cutSteps[selected]:Array.from({length:lastPosition()+1},(_,i)=>i);
 const choices=sequence.filter(i=>direction>0?i>position():i<position());
 if(choices.length)move(direction>0?choices[0]:choices.at(-1));
 if(position()===lastPosition())pause();
}
el('first').onclick=()=>{pause();move(0);};
el('last').onclick=()=>{pause();move(lastPosition());};
el('previous').onclick=()=>{pause();step(-1);};el('next').onclick=()=>{pause();step(1);};
el('slider').oninput=el('round').oninput=event=>{pause();move(Math.max(0,Math.min(lastPosition(),Number(event.target.value))));};
el('budget').onchange=()=>{pause();selected=Number(el('budget').value);current=scenarios[selected].history.length;render();};
el('mode').onchange=()=>{pause();current=budgetMode()?scenarios[selected].history.length:0;render();};
el('outer').onchange=el('before').onchange=render;el('speed').onchange=pause;
el('play').onclick=()=>{if(timer){pause();return;}if(position()===lastPosition())move(0);
 el('play').textContent='暂停';timer=setInterval(()=>step(1),Number(el('speed').value));};
render();
</script></body></html>"""
