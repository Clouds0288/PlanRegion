"""Independent AC checks; an infeasible tested design is not a global proof."""
from run_domain import ROOT, threadpool_limits
from pathlib import Path
from time import perf_counter
import json
import warnings
import numpy as np
from Network.jiangkou import Jiangkou
from vertify import ACPowerFlow
from region import sample_region, halfspaces, contains
from experiment_runtime import save, serial

LIMIT = 3600.
DIVISIONS = 100


def currents_limit(network, choice):
    catalog = network.data['line_types']
    return np.array([catalog[e['options'][k]['line_type']]['max_i_ka']
                     for e,k in zip(network.corridors,choice)]) * network.data['line_loading_limit']


def full_flow(network, power):
    oracle = ACPowerFlow(network)
    points = np.atleast_2d(power)
    ell = np.zeros((len(points), network.n))
    stable = np.zeros(len(points), bool)
    for _ in range(160):
        with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
            P,Q,v,u = oracle.state(points,ell)
            new = (P*P+Q*Q)/u
            stable = (np.max(np.abs(new-ell),axis=1)<1e-12) & (u.min(axis=1)>0)
        if stable.all():
            break
        ell = new
        if not np.isfinite(ell).all():
            break
    P,Q,v,u = oracle.state(points,ell)
    residual = np.max(np.abs(P*P+Q*Q-u*ell),axis=1)
    return P,Q,v,u,ell,stable & np.isfinite(residual) & (residual<1e-9)


def point_metrics(network, choice, power):
    c = network.design(choice)
    P,Q,v,u,ell,converged = full_flow(c,power)
    if not converged[0]:
        return dict(power_kw=power, converged=False,
                    fixed_scheme_status=int(ACPowerFlow(c).classify(power)[0]))
    current = np.sqrt(np.maximum(ell[0],0.))*c.base/(np.sqrt(3)*c.voltage_kv*1000)
    utilization = current/currents_limit(network,choice)
    ps,qs = P[0,c.roots].sum()*c.base, Q[0,c.roots].sum()*c.base
    return dict(power_kw=power,converged=True,minimum_voltage_pu=float(np.sqrt(v.min())),
                minimum_voltage_node=int(c.nodes[np.argmin(v)]),
                maximum_current_utilization=float(utilization.max()),
                most_loaded_corridor=c.corridors[int(utilization.argmax())]['id'],
                source_kw=ps,source_kvar=qs,source_kva=float(np.hypot(ps,qs)),
                loss_kw=float(np.sum(ell*c.r)*c.base),
                equality_residual=float(np.max(np.abs(P*P+Q*Q-u*ell))),
                same_contract_violation=float(ACPowerFlow(c).violation(P,Q,v).max()))


def boundary(network, choice, directions, strict=False):
    c=network.design(choice)
    oracle=ACPowerFlow(c)
    lower=np.zeros(len(directions)); upper=np.full(len(directions),network.power_limit)
    imax=currents_limit(network,choice)
    unknown=0
    for _ in range(22):
        middle=(lower+upper)/2
        status,ell=oracle.classify(directions*middle[:,None],return_currents=True)
        unknown+=int(np.sum(status==0))
        feasible=status==1
        if strict:
            current=np.sqrt(np.maximum(ell,0))*c.base/(np.sqrt(3)*c.voltage_kv*1000)
            feasible &= np.max(current/imax,axis=1)<=1+1e-9
        lower[feasible]=middle[feasible]; upper[~feasible]=middle[~feasible]
    return dict(inner_kw=directions*lower[:,None], outer_kw=directions*upper[:,None],
                radius_tolerance_kw=float(np.max(upper-lower)),design='all_max',unknown_checks=unknown)


def nodal_check(network, cases):
    import pandapower as pp
    from pandapower.converter.pypower.from_ppc import from_ppc
    output=[]
    for name,choice,power in cases:
        c=network.design(choice)
        P,Q,v,u,ell,converged=full_flow(c,power)
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',category=FutureWarning)
                net=from_ppc(c.ppc(power),f_hz=50)
            pp.runpp(net,algorithm='nr',tolerance_mva=1e-10,numba=False,max_iteration=100)
            vm=net.res_bus.loc[list(c.nodes),'vm_pu'].to_numpy()
            ps=float(net.res_ext_grid.p_mw.sum()*1000)
            qs=float(net.res_ext_grid.q_mvar.sum()*1000)
            output.append(dict(name=name,power_kw=power,converged=bool(net.converged),
                minimum_voltage_pu=float(vm.min()),source_kva=float(np.hypot(ps,qs)),
                voltage_difference_pu=float(np.max(np.abs(vm-np.sqrt(v[0])))) if converged[0] else None,
                source_p_difference_kw=float(abs(ps-P[0,c.roots].sum()*1000)) if converged[0] else None))
        except Exception as exc:
            output.append(dict(name=name,power_kw=power,converged=False,error=str(exc)))
    return dict(pandapower_version=pp.__version__,cases=output)


def main():
    started=perf_counter(); deadline=started+LIMIT
    c=Jiangkou(load_nodes=('B000025','B000078'))
    final=ROOT/'results'/'domain.json'
    domain=json.loads(final.read_text(encoding='utf-8'))
    if domain['status']=='partial_error':
        # Preliminary AC work remains valid while the global search continues.
        domain=json.loads((ROOT/'results'/'domain_checkpoint.json').read_text(encoding='utf-8'))
    bounds=np.repeat(c.power_limit,2)
    axis=(np.arange(DIVISIONS)+.5)*c.power_limit/DIVISIONS
    xx,yy=np.meshgrid(axis,axis,indexing='xy')
    points=np.c_[xx.ravel(),yy.ravel()]
    socp=sample_region(domain,points,bounds)
    max_choice=np.array([len(o.cost)-1 for o in c.line_options])
    choices=[max_choice]
    for record in domain['certificates']:
        if not record['inner']:
            continue  # Search candidates without any feasible witness are not certified designs.
        choice=np.asarray(record['choice'],int)
        if not any(np.array_equal(choice,existing) for existing in choices):
            choices.append(choice)
    ac=np.zeros(len(points),np.int8); strict=ac.copy()
    tested_negative=np.zeros(len(points),int)
    owner=np.full(len(points),-1,int)
    support_count=support_failed=support_unknown=0
    vertex_count=vertex_failed=vertex_unknown=0
    vertex_thermal_failed=0; vertex_thermal_max=0.
    for k,choice in enumerate(choices):
        oracle=ACPowerFlow(c.design(choice)); imax=currents_limit(c,choice)
        for start in range(0,len(points),512):
            if perf_counter()>=deadline:
                break
            ids=np.arange(start,min(start+512,len(points)))
            status,ell=oracle.classify(points[ids],return_currents=True)
            current=np.sqrt(np.maximum(ell,0))*c.base/(np.sqrt(3)*c.voltage_kv*1000)
            good=status==1; thermal=np.max(current/imax,axis=1)<=1+1e-9
            ac[ids[good]]=1; owner[ids[good & (owner[ids]<0)]]=k
            strict[ids[good & thermal]]=1
            tested_negative[ids[status==-1]]+=1
        for record in domain['certificates']:
            if not np.array_equal(choice,record['choice']):
                continue
            vertices=np.asarray(record['inner'])*bounds
            if not len(vertices):
                continue
            vertex_status,vertex_ell=oracle.classify(vertices,return_currents=True)
            vertex_count+=len(vertices)
            vertex_failed+=int(np.sum(vertex_status==-1));vertex_unknown+=int(np.sum(vertex_status==0))
            vertex_utilization=np.max(np.sqrt(np.maximum(vertex_ell,0))*c.base/
                                      (np.sqrt(3)*c.voltage_kv*1000)/imax,axis=1)
            vertex_thermal_failed+=int(np.sum((vertex_status==1)&(vertex_utilization>1+1e-9)))
            if np.any(vertex_status==1):
                vertex_thermal_max=max(vertex_thermal_max,float(vertex_utilization[vertex_status==1].max()))
            inside=contains(points/bounds,halfspaces(np.asarray(record['inner'])))
            ids=np.flatnonzero(inside)
            for start in range(0,len(ids),512):
                status=oracle.classify(points[ids[start:start+512]])
                support_count+=len(status)
                support_failed+=int(np.sum(status==-1));support_unknown+=int(np.sum(status==0))
        np.savez_compressed(ROOT/'results'/'ac_grid.npz',points_kw=points,axis_kw=axis,
            socp=socp,ac_verified=ac,ac_strict_current=strict,owner=owner,
            all_tested_negative=tested_negative==len(choices),choices=np.asarray(choices))
        print(f'AC {k+1}/{len(choices)} designs, {perf_counter()-started:.1f}s',flush=True)
    directions=np.c_[np.linspace(1,0,101),np.linspace(0,1,101)]
    ac_boundary=boundary(c,max_choice,directions)
    strict_boundary=boundary(c,max_choice,directions,strict=True)
    save('ac_boundary.json',dict(same_contract=ac_boundary,strict_current=strict_boundary))
    nominal=c.original_p[c.selected]
    probes=[('Nominal',max_choice,nominal),('Zero variable loads',max_choice,np.zeros(2))]
    for i in (0,25,50,75,100):
        probes.append((f'Boundary direction {i}',max_choice,ac_boundary['inner_kw'][i]))
    probes.extend([('Just outside AC boundary',max_choice,ac_boundary['outer_kw'][50]*1.002)])
    save('ac_points.json',dict(all_max=[dict(name=name,**point_metrics(c,ch,p)) for name,ch,p in probes],
                             original=point_metrics(c,np.zeros(c.n,int),nominal)))
    cross=nodal_check(c,probes)
    save('ac_nodal_crosscheck.json',cross)
    complete=perf_counter()<deadline
    summary=dict(status='sampled' if complete else 'time_limit',seconds=perf_counter()-started,
        limit_seconds=LIMIT,grid_shape=[DIVISIONS,DIVISIONS],grid_spacing_kw=c.power_limit/DIVISIONS,
        input_domain_status=domain['status'],tested_designs=len(choices),
        samples=len(points),socp_inner=int(np.sum(socp==1)),socp_unknown=int(np.sum(socp==0)),
        ac_verified=int(np.sum(ac==1)),ac_strict_current=int(np.sum(strict==1)),
        socp_inside_ac_verified=int(np.sum((socp==1)&(ac==1))),
        socp_inside_ac_not_verified=int(np.sum((socp==1)&(ac==0))),
        ac_verified_outside_inner=int(np.sum((socp!=1)&(ac==1))),
        ac_verified_outside_global_outer=int(np.sum((socp==-1)&(ac==1))),
        socp_inside_failing_strict_current=int(np.sum((socp==1)&(ac==1)&(strict==0))),
        matched_scheme_checks=support_count,matched_scheme_failures=support_failed,
        matched_scheme_unknown=support_unknown,
        inner_vertices_checked=vertex_count,inner_vertices_failed=vertex_failed,inner_vertices_unknown=vertex_unknown,
        inner_vertices_exceeding_current_threshold=vertex_thermal_failed,
        inner_vertex_maximum_current_utilization=vertex_thermal_max,
        caveat='AC feasible means at least one tested design passes. Other samples are not globally AC-infeasible proofs. Ray curves refer to the all-max design.')
    save('ac_summary.json',summary)
    print(json.dumps(serial(summary),ensure_ascii=False),flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):
        main()
