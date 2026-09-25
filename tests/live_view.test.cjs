// No browser or solver: exercise the actual viewer script and numeric display helpers.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');

class Element {
  constructor() { this.value = ''; this.checked = true; this.children = []; this.attrs = {}; this.scrollTop = 0; this.selectedOptions = [{text: ''}]; this.classList = {toggle() {}}; }
  append(...items) { this.children.push(...items); }
  replaceChildren(...items) { this.children = items; }
  add(item) { this.children.push(item); }
  setAttribute(key, value) { this.attrs[key] = value; }
  on(name, handler) { (this.handlers ||= {})[name] = handler; }
}
const elements = new Map(), get = id => { if (!elements.has(id)) elements.set(id, new Element()); return elements.get(id); };
for (const id of ['method', 'budget']) get(id).value = 'auto';
get('cut-stage').value = 'cut';
const document = {getElementById: get, createElement: () => new Element()};
const Plotly = {
  async react(id, data, layout) { Object.assign(get(id), {data, layout}); },
  async relayout(id, update) { get(id).layout.scene.camera = update['scene.camera']; await get(id).handlers?.plotly_relayout(update); }
};
const context = vm.createContext({document, Plotly, window: {Plotly}, location: {protocol: 'file:'}, Option: function(text, value) { return {text, value}; }, setTimeout, clearTimeout, console});
const html = fs.readFileSync(path.join(__dirname, '..', 'live_view.html'), 'utf8');
for (const script of html.matchAll(/<script>([\s\S]*?)<\/script>/g)) vm.runInContext(script[1], context);
const run = code => vm.runInContext(code, context);
const plain = value => JSON.parse(JSON.stringify(value));

(async () => {
  // Auxiliary box/plane cannot dwarf the actual region; each axis fits its own data.
  run('var fitted=viewAxes([boxWire([150,150,150]),{name:"region",x:[0,10],y:[0,20],z:[0,30]},{name:"选中方案条件切面",x:[150],y:[150],z:[150]}],[1,2,3],[150,150,150])');
  assert.deepEqual(plain(run('[fitted.xaxis.range,fitted.yaxis.range,fitted.zaxis.range]')), [[0,10.6],[0,21.2],[0,31.8]]);
  assert.deepEqual(plain(run('viewAxes([], [1,2,3], [150,120,100]).yaxis.range')), [0,120]);
  assert.deepEqual(plain(run('viewAxes([{x:[0],y:[0],z:[0]}], [1,2,3], [150,150,150]).xaxis.range')), [0,.06]);
  // Coefficients must be conditioned on the selected plan, not the generating plan.
  run('var cut={joint_coefficients:[-3,1,0,0,2],constant:-1,normal:[1,0,0],choice:{e:"a"},number:1}');
  assert.equal(run('conditionCut(cut,{x:[0],choice:{e:"b"}}).constant'), -3);
  assert.equal(run('conditionCut(cut,{x:[1],choice:{e:"a"}}).constant'), -1);
  assert.equal(run('conditionCut(cut,{choice:{e:"b"}})'), null);
  assert.equal(run('conditionCut(cut,{choice:{e:"a"}}).constant'), -1);
  // Clipping uses kW correctly even when the three axes have different scales.
  run('var bounds=[2,3,5],cube=boxMesh(bounds),sliced=sliceMesh(cube,{constant:1,normal:[-1,0,0]},bounds)');
  assert(run('sliced.kept.vertices.every(p=>p[0]<=1+1e-10)'));
  assert(run('sliced.plane.vertices.every(p=>Math.abs(p[0]-1)<1e-10)'));
  assert.equal(run('sliceMesh(cube,{constant:-1,normal:[0,0,0]},bounds).kept.vertices.length'), 0);
  assert.equal(run('sliceMesh(cube,{constant:1,normal:[0,0,0]},bounds).kept.vertices.length'), 8);
  assert.equal(run('sliceMesh(cube,{constant:-20,normal:[1,1,1]},bounds).kept.vertices.length'), 0);
  assert(run('coveredByOne(boxMesh([1,1,1]),cube,bounds)'));
  assert(!run('coveredByOne(cube,boxMesh([1,1,1]),bounds)'));
  assert(!run('coveredByOne(cube,{vertices:[[0,0,0]],faces:[]},bounds)'));
  // Empty final-result placeholders and genuine empty candidate domains are distinct.
  run('var empty={vertices:[],faces:[]},g={choice:{e:"a"},inner:empty,outer:empty}');
  assert.equal(run('schemeStatus(g,false)[0]'), 'empty');
  assert.equal(run('schemeStatus(g,true)[0]'), 'unknown');
  // An overlapping union covering the vertices is not a single-convex-set certificate.
  assert(!run('coveredByOne(cube,boxMesh([1.2,3,5]),bounds)'));
  assert(!run('coveredByOne(cube,{vertices:boxMesh([1.2,3,5]).vertices.map(p=>[p[0]+.8,p[1],p[2]]),faces:cube.faces},bounds)'));
  // A lower-dimensional surviving set is not classified as empty.
  run('g.outer={vertices:[[0,0,0]],faces:[]}');
  assert.equal(run('schemeStatus(g,false)[0]'), 'unknown');

  run(`var first={x:[1],choice:{e:'a'},cost:0,inner:empty,outer:cube};
       var second={x:[0],choice:{e:'b'},cost:0,inner:empty,outer:cube};
       var frame0={event:'point',method:'linear',phase:'linear',budget_index:0,budgets:[0],bounds,load_nodes:[1,2,3],geometry:[first,second],choice:first.choice};
       var frame1={event:'cut',latest_cut:cut,geometry:[{...first,outer:sliceMesh(cube,{constant:-1,normal:[1,0,0]},bounds).kept},{...second,outer:empty}]};
       ingest({...frame0,...frame1,history:[{id:0,elapsed:0,patch:frame0},{id:1,elapsed:1,patch:frame1}]},true);`);
  await run('seek(1)');
  assert.equal(run('cutHistory(state).cut.number'), 1);
  assert.deepEqual(plain(run('[...cutChanges(cutHistory(state)).values()]')), ['收紧', '切空']);
  assert.equal(get('scheme-list').children.length, 2);
  await run('selectScheme(choiceKey(second.choice))');
  assert(get('scheme-status').textContent.includes('已证不可行'));
  assert(get('cut-equation').textContent.includes('-3.0000'));
  const actualCamera = {eye:{x:2,y:-1,z:.8},up:{x:0,y:0,z:1},center:{x:0,y:0,z:0}};
  await Plotly.relayout('chart', {'scene.camera':actualCamera});
  assert.deepEqual(plain(get('detail-chart').layout.scene.camera), actualCamera);
  await run('seek(0)');
  assert.deepEqual(plain(get('chart').layout.scene.camera), actualCamera);
  assert(!get('scheme-status').textContent.includes('已证不可行'));
  await run('seek(1)');
  get('sync-camera').checked = false;
  await Plotly.relayout('chart', {'scene.camera':{...actualCamera,eye:{x:3,y:2,z:1}}});
  assert.deepEqual(plain(get('detail-chart').layout.scene.camera), actualCamera);
  await get('reset-camera').onclick();
  assert.deepEqual(plain(get('chart').layout.scene.camera), plain(run('defaultCamera')));
  run("ingest({...frame0,phase:'socp',history:[{id:2,elapsed:2,patch:{event:'phase_start',phase:'socp',latest_cut:null,geometry:[]}}]})");
  await run('seek(2)');
  assert.equal(run('cutHistory(state)'), null);
  assert.equal(get('scheme-list').children.length, 0);
  console.log('Viewer geometry, evidence states, replay, selection and camera tests passed.');
})().catch(error => { console.error(error); process.exitCode = 1; });
