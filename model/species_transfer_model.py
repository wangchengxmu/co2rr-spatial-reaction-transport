"""Effective-shear sensitivity: species-dependent passive D/k, not a thin slab.

Original chemistry and previously saved calculations are not modified.
See ../docs/MODEL.md for assumptions and validity limits.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import argparse, json, hashlib, time, traceback
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.linalg import solve_banded
import corrected_numerics as n

m=n.m
HERE=Path(__file__).resolve().parent
ROOT=HERE/'species_transfer_screen'

@dataclass(frozen=True)
class Config(m.Config):
    delta_CO2_um:float=10.
    n_catalyst:int=48
    n_external:int=96
    n_axial:int=24
    allocation:str='CO2_weighted'

def hydro(cfg):
    length=np.sqrt(cfg.area_m2)
    gamma=m.D[0]*length/(.807*cfg.delta_CO2_um*1e-6)**3
    km=.807*(m.D**2*gamma/length)**(1/3)
    real_U=cfg.flow_mL_min*1e-6/60/(length*cfg.channel_height_m)
    return {'length_m':float(length),'effective_shear_s-1':float(gamma),
      'mean_velocity_m_s':float(gamma*cfg.channel_height_m/6),
      'actual_flow_reference_shear_s-1':float(6*real_U/cfg.channel_height_m),
      'passive_delta_CO2_m':cfg.delta_CO2_um*1e-6,'passive_k_CO2_m_s':float(km[0]),
      'k_species_m_s':dict(zip(m.SP,km.tolist())),
      'delta_species_m':dict(zip(m.SP,(m.D/km).tolist())),
      'definition':'Effective linear shear sensitivity; not pump-derived velocity. Passive species D/k lengths, common extended liquid domain.'}

m.hydrodynamics=hydro
# Reuse audited finite-volume chemistry, changing only the velocity profile.
import inspect
problem_source=inspect.getsource(n.station_problem)
old="u=6*h['mean_velocity_m_s']*(y/cfg.channel_height_m)*(1-y/cfg.channel_height_m)"
assert problem_source.count(old)==1
exec(problem_source.replace(old,"u=h['effective_shear_s-1']*y"),n.__dict__)

def jacobian(fun,t,v,ncat):
    """Local coloured differences plus exact nonlocal catalyst-CO2 columns.

    Normalized CO2 weights and cumulative ionic current couple these columns
    to all cells. A purely tridiagonal Jacobian is wrong for this model.
    """
    count=len(v)//6;f0=fun(t,v);rr=[];cc=[];vv=[]
    for b in range(6):
        for offset in range(3):
            cells=np.arange(offset,count,3)
            if b==0:cells=cells[cells>=ncat]
            if not len(cells):continue
            cols=b*count+cells;steps=1e-6*np.maximum(abs(v[cols]),1.)
            trial=v.copy();trial[cols]+=steps;change=fun(t,trial)-f0
            for cell,col,step in zip(cells,cols,steps):
                rows=(np.arange(6)[:,None]*count+np.arange(max(0,cell-1),min(count,cell+2))).ravel()
                rr.extend(rows);cc.extend([col]*len(rows));vv.extend(change[rows]/step)
    for col in range(ncat):
        step=1e-6*max(abs(v[col]),1.);trial=v.copy();trial[col]+=step
        change=(fun(t,trial)-f0)/step;rows=np.flatnonzero(change)
        rr.extend(rows);cc.extend([col]*len(rows));vv.extend(change[rows])
    return coo_matrix((vv,(rr,cc)),shape=(6*count,6*count)).tocsc()

original_ivp=n.solve_ivp
ACTIVE_NCAT=48
def bounded_ivp(*args,**kwargs):
    kwargs.pop('jac_sparsity',None)
    kwargs['jac']=lambda t,v:jacobian(args[0],t,v,ACTIVE_NCAT)
    return original_ivp(*args,**kwargs)
n.solve_ivp=bounded_ivp
m.solve_station=n.solve_station

def passive(cfg,diffusivity=None,nx=None,ny=None,outer=None):
    """Independent scalar FV check, same linear flow, no reaction/current."""
    from dataclasses import replace
    c=replace(cfg,n_external=ny or cfg.n_external,outer_domain_m=outer or cfg.outer_domain_m)
    edges,z,w,_=m.mesh(c);edges=edges[c.n_catalyst:]-c.catalyst_length_m
    z=z[c.n_catalyst:]-c.catalyst_length_m;w=w[c.n_catalyst:]
    D=float(m.D[0] if diffusivity is None else diffusivity)
    hy=hydro(c);u=hy['effective_shear_s-1']*z;L=hy['length_m']
    xe=L*np.linspace(0,1,(nx or c.n_axial)+1)**1.3
    g=D/np.diff(z);g0=D/z[0];gn=D/(edges[-1]-z[-1]);state=np.ones(len(z));flux=0.
    for dx in np.diff(xe):
        left=np.r_[g0,g]/w;right=np.r_[g,gn]/w
        diag=left+right+u/dx;ab=np.zeros((3,len(z)))
        ab[0,1:]=-right[:-1];ab[1]=diag;ab[2,:-1]=-left[1:]
        rhs=u/dx*state;rhs[-1]+=right[-1]
        state=solve_banded((1,1),ab,rhs);flux+=g0*state[0]*dx/L
    expected=.807*(D*D*hy['effective_shear_s-1']/L)**(1/3)
    return dict(delta_CO2_target_um=cfg.delta_CO2_um,D_m2_s=D,nx=len(xe)-1,ny=len(z),
      outer_um=c.outer_domain_m*1e6,k_m_s=float(flux),k_Leveque_m_s=float(expected),
      delta_FV_um=float(D/flux*1e6),relative_k_error=float(flux/expected-1))

def hashes():
    paths=[Path(__file__),HERE/'corrected_numerics.py',HERE/'current_constrained_model.py',HERE.parent/'fe_audit/faradaic_inputs.csv']
    return {p.relative_to(HERE.parent).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2),encoding='utf-8');tmp.replace(path)

def run(cfg,out):
    global ACTIVE_NCAT
    ACTIVE_NCAT=cfg.n_catalyst;out.mkdir(parents=True,exist_ok=True)
    source_hashes=hashes();write(out/'config.json',asdict(cfg))
    stations=[];_,z,_,_=m.mesh(cfg);L=hydro(cfg)['length_m']
    xe=L*np.linspace(0,1,cfg.n_axial+1)**1.3
    upstream=np.repeat(m.bulk(cfg)[:,None],len(z),axis=1)
    for ix,dx in enumerate(np.diff(xe)):
        s=n.solve_station(cfg,upstream,dx);s['x_m']=(xe[ix]+xe[ix+1])/2
        if s['components'][0].min()<0:raise RuntimeError('Negative raw CO2')
        stations.append(s);upstream=s['components']
        write(out/'progress.json',{'finished_station':ix+1,'total_stations':cfg.n_axial,'diagnostics':s['diagnostics']})
        print(json.dumps({'station':ix+1,'total':cfg.n_axial,**s['diagnostics']}),flush=True)
    raw=out/'raw.npz'
    np.savez_compressed(raw,components=np.array([s['components'] for s in stations]),
        upstream=np.array([s['upstream'] for s in stations]),z_m=z,dx_m=np.diff(xe),x_m=[s['x_m'] for s in stations])
    gas=m.isotope_channel(cfg,stations,1.,0.);liquid=m.isotope_channel(cfg,stations,0.,1.)
    closure=max(float(abs(g['components']+l['components']-np.vstack([s['c'][0],s['c'][1]+s['c'][2]])).max()) for s,g,l in zip(stations,gas,liquid))
    tracer_min=min(float(t['components'].min()) for t in gas+liquid)
    fraction_error=max(float(np.max(np.maximum(-t['fraction_CO2'],t['fraction_CO2']-1))) for t in gas+liquid)
    assert closure<1e-5 and tracer_min>=-1e-8 and fraction_error<1e-6,(closure,tracer_min,fraction_error)
    # Recompute every station from saved raw arrays before accepting outputs.
    reload=np.load(raw)
    for ix,s in enumerate(stations):
        ok,_=n.station_problem(cfg,reload['upstream'][ix],float(reload['dx_m'][ix]))[1](reload['components'][ix])
        assert ok
    m.OUT=out;report=m.save(cfg,stations,'case')
    frames=[];r=m.rates(cfg)
    for s,g,l in zip(stations,gas,liquid):
        frame=pd.DataFrame({'z_nm':z*1e9,'region':np.where(z<cfg.catalyst_length_m,'catalyst','liquid'),
          'x_mm':s['x_m']*1e3,'dx_m':s['dx'],
          'carbon_current_mA_cm2_per_nm':r['carbon_electrons']*s['weights']*m.FARADAY*1e-10,
          'total_current_mA_cm2_per_nm':s['electron_source']*m.FARADAY*1e-10,
          'gas_CO2_mM':g['components'][0],'liquid_CO2_mM':l['components'][0],
          'gas_ionic_DIC_mM':g['components'][1],'liquid_ionic_DIC_mM':l['components'][1]})
        frames.append(frame)
    pd.concat(frames).to_csv(out/'current_and_tracer_axial.csv',index=False)
    df=pd.read_csv(out/'case_mean_profile.csv')
    for col in ['carbon_current_mA_cm2_per_nm','total_current_mA_cm2_per_nm']:
        df[col]=sum(f[col].to_numpy()*s['dx']/L for f,s in zip(frames,stations))
    df['noncarbon_current_mA_cm2_per_nm']=df.total_current_mA_cm2_per_nm-df.carbon_current_mA_cm2_per_nm
    df.to_csv(out/'case_mean_profile.csv',index=False)
    report.update(transport_model='extended_liquid_effective_linear_shear',
      all_raw_restart_checks_passed=True,max_pointwise_isotope_closure_mM=closure,
      min_isotope_component_mM=tracer_min,max_isotope_fraction_violation=fraction_error,
      raw_CO2_min_mM=min(float(s['components'][0].min()) for s in stations),
      passive_benchmark=passive(cfg),input_hashes=source_hashes)
    assert source_hashes==hashes()
    write(out/'validated.json',report)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--delta',type=float,default=10);p.add_argument('--chi',type=float,default=.25)
    p.add_argument('--nc',type=int,default=48);p.add_argument('--ny',type=int,default=96);p.add_argument('--nx',type=int,default=24)
    p.add_argument('--outer',type=float,default=500);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();cfg=Config(delta_CO2_um=args.delta,chi=args.chi,n_catalyst=args.nc,n_external=args.ny,n_axial=args.nx,outer_domain_m=args.outer*1e-6)
    t0=time.monotonic()
    try:
        report=run(cfg,args.out);result={'status':'passed','summary':report['summary']}
    except Exception as exc:
        result={'status':'not_accepted','reason':str(exc),'traceback':traceback.format_exc()}
    result['elapsed_s']=time.monotonic()-t0;write(args.out/'result.json',result)
    print(json.dumps(result),flush=True)
