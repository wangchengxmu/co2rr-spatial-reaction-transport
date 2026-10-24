"""Numerical corrections only; original model and old results remain immutable.

Import a private copy of the original model, tightening its algebraic closure.
Direct steady solves supplement time relaxation; acceptance never uses solver
success alone. Raw states are the canonical restart, not species CSVs.
"""
from pathlib import Path
import importlib.util,sys
import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('_corrected_base',HERE/'current_constrained_model.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
# Exact, recorded numerical change, no chemical equations changed.
source=(HERE/'current_constrained_model.py').read_text(encoding='utf-8')
closure=source[source.index('def speciate('):source.index('def bulk(')]
assert closure.count('2e-10')==2
exec(closure.replace('2e-10','2e-12'),m.__dict__)
TOL=2e-6

def station_problem(cfg,upstream,dx):
    e,z,w,d=m.mesh(cfg);n=len(z);h=m.hydrodynamics(cfg)
    y=np.maximum(z-cfg.catalyst_length_m,0)
    u=6*h['mean_velocity_m_s']*(y/cfg.channel_height_m)*(1-y/cfg.channel_height_m)
    eps=np.where(z<cfg.catalyst_length_m,cfg.porosity,1.)
    k=m.constants(cfg);r=m.rates(cfg);kh=cfg.hydration_s;ko=cfg.hydroxylation_M_s/1000
    def rhs(t,flat):
        v=flat.reshape(6,n);c=m.speciate(v,cfg)
        f,_,weights,se=m.fluxes(c,cfg,e,z,w,d)
        dv=-(m.W@(f[:,1:]-f[:,:-1]))/w-u/dx*(v-upstream)
        net=eps*((kh+ko*c[4])*c[0]-(kh/k['K1']*c[3]+ko*k['KW']/k['K1'])*c[1])
        dv[0]-=net+r['carbon']*weights;dv[1]+=net;dv[2]+=se
        dv[4]+=r['formate']*weights;dv[5]+=r['acetate']*weights
        return (dv/eps).ravel()
    def assess(v):
        c=m.speciate(v,cfg);der=rhs(0,v.ravel()).reshape(v.shape)
        f,psi,weights,se=m.fluxes(c,cfg,e,z,w,d)
        carbon=f[0,0]-sum(f[:3,-1])-r['carbon']-np.dot(u/dx,((v[0]+v[1])-(upstream[0]+upstream[1]))*w)
        alkflux=m.W[2]@f
        alk=alkflux[0]-alkflux[-1]+r['electron']-np.dot(u/dx,(v[2]-upstream[2])*w)
        diag={'steady_scaled':float(np.max(abs(der)/np.maximum(abs(v),1))),
              'carbon_balance_mol_m2_s':float(carbon),'alkalinity_balance_mol_m2_s':float(alk),
              'max_charge_mM':float(abs(m.Z@c).max()),'max_water_relative_error':float(abs(c[3]*c[4]/k['KW']-1).max()),
              'max_ionic_drift_step':float(abs(m.Z[:,None]*psi[None,1:]).max()),
              'min_component_mM':float(v[[0,1,3,4,5]].min()),'gas_C_flux_mol_m2_s':float(f[0,0]),
              'liquid_edge_DIC_flux_mol_m2_s':float(sum(f[:3,cfg.n_catalyst])), 'max_CO2_mM':float(c[0].max()),
              'carbon_source_integral':float(weights@w),'electron_source_integral':float(se@w),
              'max_cell_alkalinity_balance_mol_m2_s':float(abs(der[2]*eps*w).max())}
        passed=(np.isfinite(v).all() and np.isfinite(der).all() and diag['steady_scaled']<TOL
                and diag['min_component_mM']>=-1e-8 and c.min()>=-1e-8
                and abs(carbon)<1e-10 and abs(alk)<1e-10 and diag['max_charge_mM']<1e-8
                and diag['max_water_relative_error']<1e-10 and diag['max_ionic_drift_step']<2
                and abs(weights@w-1)<1e-10 and abs(se@w-r['electron'])<1e-12)
        return passed,dict(components=v.copy(),c=c,flux=f,psi=psi,weights=weights,electron_source=se,
                          z=z,w=w,d=d,edges=e,u=u,dx=dx,diagnostics=diag,upstream=upstream.copy())
    return rhs,assess

def solve_station(cfg,upstream,dx,initial=None):
    rhs,assess=station_problem(cfg,upstream,dx)
    state=(upstream if initial is None else initial).copy();shape=state.shape;attempts=[]
    def polish(v):
        scale=np.maximum(abs(v),1).ravel()
        try:
            sol=root(lambda x:rhs(0,x*scale)/scale,v.ravel()/scale,method='hybr',options={'xtol':1e-10,'eps':1e-12,'maxfev':12000,'factor':.1})
        except (RuntimeError,ValueError,FloatingPointError) as exc:
            attempts.append({'method':'direct_root','rejected_trial':str(exc)})
            return v.copy()
        out=(sol.x*scale).reshape(shape)
        attempts.append({'method':'direct_root','reported_success':bool(sol.success),'nfev':int(sol.nfev),
                         'max_correction_mM':float(abs(out-v).max())})
        return out
    for phase in ['initial','root','relax']:
        if phase=='root':
            candidate=polish(state)
            ok,result=assess(candidate)
            if ok:result['numerical_attempts']=attempts;return result
            # Do not integrate an invalid root iterate.
            if candidate[[0,1,3,4,5]].min()>=-1e-8 and np.isfinite(candidate).all():state=candidate
        if phase=='relax':
            for stop in [.02,.2,2.,20.,100.]:
                def depletion(t,v):return np.min(v[:shape[1]])+1e-8
                depletion.terminal=True;depletion.direction=-1
                sol=solve_ivp(rhs,(0,stop),state.ravel(),method='BDF',rtol=1e-9,atol=2e-10,
                              jac_sparsity=m.sparsity(shape[1],cfg.n_catalyst),events=depletion)
                if not sol.success or sol.status==1:raise RuntimeError('Relaxation failed or raw CO2 depletion: '+str(sol.message))
                state=sol.y[:,-1].reshape(shape);ok,result=assess(state)
                attempts.append({'method':'relax','duration_s':stop,'residual':result['diagnostics']['steady_scaled']})
                if not ok and stop>=2:
                    candidate=polish(state);ok,checked=assess(candidate)
                    if ok:result=checked
                if ok:result['numerical_attempts']=attempts;return result
        else:
            ok,result=assess(state)
            if ok:result['numerical_attempts']=attempts;return result
    raise RuntimeError('Station failed unchanged steady-state gate: '+str(result['diagnostics']))
