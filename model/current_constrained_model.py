"""Prescribed-current, FE-resolved porous-catalyst / laminar-channel model.

No electrochemical rate constant is fitted. The integrated carbon sink is
fixed by j*sum(FE*nC/ne)/F at every axial station. Its transverse distribution
is an explicitly assumed CO2-concentration weighting, not measured kinetics.
Fast proton-transfer constraints do NOT impose CO2/DIC equilibrium.

Coordinates: z=0 gas contact; z=Lc catalyst/catholyte boundary; y=z-Lc in
the liquid. x follows the channel. A steady upstream march resolves axial
convection and transverse diffusion/migration/reaction. The outer domain is
a numerical truncation, not a stagnant-film thickness.
Units: mol/m3 (numerically mM), m, s. Ideal activities remain an assumption.
"""
from dataclasses import dataclass, asdict, replace
from pathlib import Path
import argparse, json, hashlib, time
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
FARADAY=96485.33212
SP=['CO2','HCO3','CO3','H','OH','Cs','HSO4','SO4','HCOOH','HCOO','CH3COOH','CH3COO']
Z=np.array([0,-1,-2,1,-1,1,-1,-2,0,-1,0,-1.],float)
# Dilute-aqueous transport proxies. See parameter provenance registry.
D=np.array([2.04527,1.10059,.804197,9.2983996,5.3054587,2.07,1.33,1.07,1.5,1.454,1.21,1.089])*1e-9
# Conserved/component concentrations: CO2, bicarbonate+carbonate, total
# alkalinity (including organic anions), total sulfate, total formate, total acetate.
W=np.zeros((6,len(SP)))
W[0,0]=1;W[1,1:3]=1
W[2,[1,2,3,4,7,9,11]]=[1,2,-1,1,1,1,1]
W[3,6:8]=1;W[4,8:10]=1;W[5,10:12]=1

@dataclass(frozen=True)
class Config:
    current_mA_cm2:float=300.
    catalyst_length_m:float=200e-9
    chi:float=.25
    porosity:float=.5
    flow_mL_min:float=10.
    area_m2:float=.8e-4
    channel_height_m:float=.002
    outer_domain_m:float=.0005
    n_catalyst:int=24
    n_external:int=56
    n_axial:int=10
    pH_bulk:float=3.4
    Cs_bulk_M:float=1.
    CO2_star_mM:float=33.
    pKa1:float=6.35
    pKa2:float=10.33
    pKaS:float=1.99
    pKa_formic:float=3.75
    pKa_acetic:float=4.756
    pKw:float=14.
    hydration_s:float=.037
    hydroxylation_M_s:float=12100.
    unresolved_FE_mode:str='noncarbon'
    allocation:str='CO2_weighted'

@lru_cache(maxsize=1)
def fe_inputs():
    path=HERE.parent/'fe_audit/faradaic_inputs.csv'
    if not path.exists():raise FileNotFoundError(path)
    return pd.read_csv(path)

def rates(cfg):
    """Mean ordinary-CO2RR FE prescription, not isotope-run-specific fitting."""
    # Source means recomputed from the three experimental replicates; no FE normalization.
    products=[(row['product'],float(row['FE_mean_percent']),int(row['n_carbon']),int(row['n_electron']))
              for row in fe_inputs().to_dict('records')]
    lookup={name:fe/100/ne for name,fe,nc,ne in products}
    q=cfg.current_mA_cm2*10/FARADAY
    frac=sum(x[1] for x in products)/100
    result={'electron':q,'carbon':q*sum(fe/100*nc/ne for _,fe,nc,ne in products),
      'carbon_electrons':q*sum(fe/100 for _,fe,nc,ne in products if nc),
      'formate':q*lookup['formate'],'acetate':q*lookup['acetate'],
      'FE_sum':frac,'unclosed_FE':1-frac,'products':products}
    if cfg.unresolved_FE_mode in ['CO','formate']:
        result['carbon']+=q*(1-frac)/2
        result['carbon_electrons']+=q*(1-frac)
        if cfg.unresolved_FE_mode=='formate':result['formate']+=q*(1-frac)/2
    elif cfg.unresolved_FE_mode!='noncarbon':raise ValueError(cfg.unresolved_FE_mode)
    result['noncarbon_electrons']=q-result['carbon_electrons']
    return result

def constants(cfg):
    return dict(K1=1000*10**(-cfg.pKa1),K2=1000*10**(-cfg.pKa2),
      KS=1000*10**(-cfg.pKaS),KF=1000*10**(-cfg.pKa_formic),
      KA=1000*10**(-cfg.pKa_acetic),KW=1e6*10**(-cfg.pKw))

def speciate(v,cfg):
    """Water/proton-transfer equilibrium at independent CO2 and ionic carbon.

    Electroneutrality implies Cs = alkalinity + total sulfate. The alkalinity
    definition includes SO4 relative to HSO4 and organic carboxylate anions.
    A safeguarded vector Newton solve locates the unique positive H+ root.
    """
    C,B,A,S,P,Q=v
    C=np.maximum(C,0);B=np.maximum(B,0);S=np.maximum(S,0)
    P=np.maximum(P,0);Q=np.maximum(Q,0)
    k=constants(cfg)
    lo=np.full_like(C,-40.);hi=np.full_like(C,15.)
    excess=A-S-2*B-P-Q
    guess=np.where(excess>0,np.log(k['KW']/np.maximum(excess,1e-30)),
                   np.log(np.maximum(k['KS']*(S/np.maximum(A-B,1e-12)-1),1e-12)))
    t=np.clip(guess,lo,hi)
    for _ in range(48):
        h=np.exp(t);o=k['KW']/h
        c2=B*k['K2']/(h+k['K2']);s2=S*k['KS']/(h+k['KS'])
        p=P*k['KF']/(h+k['KF']);a=Q*k['KA']/(h+k['KA'])
        res=B+c2+o-h+s2+p+a-A
        der=-c2*h/(h+k['K2'])-o-h-s2*h/(h+k['KS'])-p*h/(h+k['KF'])-a*h/(h+k['KA'])
        done=np.abs(res)<2e-10
        lo=np.where((res>0)&~done,t,lo);hi=np.where((res<0)&~done,t,hi)
        new=t-res/np.minimum(der,-1e-300)
        new=np.where((new<=lo)|(new>=hi)|~np.isfinite(new),.5*(lo+hi),new)
        if np.all(done):break
        t=np.where(done,t,new)
    else:
        # BDF trial states can cross sharp buffer transitions. Safeguarded
        # Newton may then advance too slowly; finish with bracket bisection.
        for _ in range(96):
            t=np.where(done,t,.5*(lo+hi))
            h=np.exp(t);o=k['KW']/h
            c2=B*k['K2']/(h+k['K2']);s2=S*k['KS']/(h+k['KS'])
            p=P*k['KF']/(h+k['KF']);a=Q*k['KA']/(h+k['KA'])
            res=B+c2+o-h+s2+p+a-A
            done=np.abs(res)<2e-10
            if np.all(done):break
            lo=np.where((res>0)&~done,t,lo);hi=np.where((res<0)&~done,t,hi)
        if np.max(abs(res))>1e-6:raise RuntimeError(f'proton-transfer closure {np.max(abs(res))}')
    return np.vstack([C,B-c2,c2,h,o,A+S,S-s2,s2,P-p,p,Q-a,a])

def bulk(cfg):
    k=constants(cfg);h=1000*10**(-cfg.pH_bulk);o=k['KW']/h
    c0=cfg.CO2_star_mM;c1=c0*k['K1']/h;c2=c1*k['K2']/h
    cs=cfg.Cs_bulk_M*1000
    frac=k['KS']/(h+k['KS'])
    s=(cs+h-o-c1-2*c2)/(1+frac)
    return np.array([c0,c1+c2,cs-s,s,0.,0.])

def mesh(cfg):
    left=np.linspace(0,cfg.catalyst_length_m,cfg.n_catalyst+1)
    # Resolve the near-wall chemistry without assigning D/k to a physical slab.
    t=np.linspace(0,1,cfg.n_external+1)
    y=cfg.outer_domain_m*np.expm1(5*t)/np.expm1(5)
    e=np.r_[left,cfg.catalyst_length_m+y[1:]]
    z=.5*(e[1:]+e[:-1]);w=np.diff(e)
    d=D[:,None]*np.where(z<cfg.catalyst_length_m,cfg.chi,1.)
    return e,z,w,d

def hydrodynamics(cfg):
    length=np.sqrt(cfg.area_m2);width=length;h=cfg.channel_height_m
    U=cfg.flow_mL_min*1e-6/60/(width*h);shear=6*U/h
    km=.807*(D**2*shear/length)**(1/3)
    dh=2*width*h/(width+h)
    return {'length_m':length,'width_m':width,'height_m':h,'mean_velocity_m_s':U,
      'wall_shear_s-1':shear,'Re_Dh':997*U*dh/.00089,
      'passive_k_CO2_m_s':km[0],'passive_delta_CO2_m':D[0]/km[0],
      'k_species_m_s':dict(zip(SP,km.tolist())),
      'definition':'0.807 length-average Leveque coefficient; passive smooth-wall approximation, not a reactive-film geometry'}

def spatial_sources(c,cfg,widths):
    r=rates(cfg);n=len(widths);nc=cfg.n_catalyst
    weights=np.zeros(n)
    if cfg.allocation=='uniform':weights[:nc]=1/cfg.catalyst_length_m
    else:
        integral=float(np.dot(c[0,:nc],widths[:nc]))
        if integral<=1e-30:raise RuntimeError('No CO2 for prescribed current')
        weights[:nc]=c[0,:nc]/integral
    se=r['carbon_electrons']*weights
    se[:nc]+=r['noncarbon_electrons']/cfg.catalyst_length_m
    return weights,se

def fluxes(c,cfg,e,z,w,d):
    """Central Nernst-Planck discretization; drift-step validity is audited."""
    n=len(z);nf=n+1
    f=np.zeros((len(SP),nf));psi=np.zeros(nf)
    weights,se=spatial_sources(c,cfg,w)
    charge=-np.r_[0,np.cumsum(se*w)]
    f[0,0]=d[0,0]/(.5*w[0])*(cfg.CO2_star_mM-c[0,0])
    g=1/(.5*w[:-1]/d[:,:-1]+.5*w[1:]/d[:,1:])
    cb=speciate(bulk(cfg)[:,None],cfg)[:,0]
    g=np.c_[g,d[:,-1]/(.5*w[-1])]
    cr=np.c_[c[:,1:],cb]
    cl=c
    avg=np.maximum(.5*(cr+cl),1e-25)
    fd=-g*(cr-cl)
    dp=(np.sum(Z[:,None]*fd,axis=0)-charge[1:])/np.maximum(np.sum(Z[:,None]**2*g*avg,axis=0),1e-25)
    f[:,1:]=fd-g*Z[:,None]*avg*dp
    psi[1:]=dp
    return f,psi,weights,se

def sparsity(n,nc):
    p=lil_matrix((6*n,6*n),dtype=int)
    for a in range(6):
      for b in range(6):
        for i in range(n):
          for j in range(max(0,i-1),min(n,i+2)):p[a*n+i,b*n+j]=1
        # CO2-weighted source and its cumulative ionic current have nonlocal
        # derivative coupling to all catalyst cells.
        p[a*n:a*n+n,b*n:b*n+nc]=1
    return p.tocsr()

def solve_station(cfg,upstream,dx,initial=None):
    e,z,w,d=mesh(cfg);n=len(z);hyd=hydrodynamics(cfg)
    y=np.maximum(z-cfg.catalyst_length_m,0)
    u=6*hyd['mean_velocity_m_s']*(y/cfg.channel_height_m)*(1-y/cfg.channel_height_m)
    k=constants(cfg);r=rates(cfg);kh=cfg.hydration_s;ko=cfg.hydroxylation_M_s/1000
    liquid_fraction=np.where(z<cfg.catalyst_length_m,cfg.porosity,1.)
    def rhs(t,flat):
        v=flat.reshape(6,n);c=speciate(v,cfg)
        f,_,weights,se=fluxes(c,cfg,e,z,w,d)
        dv=-(W@(f[:,1:]-f[:,:-1]))/w
        dv-=u/dx*(v-upstream)
        net=liquid_fraction*((kh+ko*c[4])*c[0]-(kh/k['K1']*c[3]+ko*k['KW']/k['K1'])*c[1])
        dv[0]-=net+r['carbon']*weights;dv[1]+=net
        dv[2]+=se;dv[4]+=r['formate']*weights;dv[5]+=r['acetate']*weights
        return (dv/liquid_fraction).ravel()
    initial=upstream if initial is None else initial
    state=initial.ravel().copy();t0=0.;steady=np.inf
    for tend in [.002,.02,.2,2.,20.,100.]:
        sol=solve_ivp(rhs,(t0,tend),state,method='BDF',rtol=1e-7,atol=2e-8,jac_sparsity=sparsity(n,cfg.n_catalyst))
        if not sol.success:raise RuntimeError(sol.message)
        state=sol.y[:,-1];t0=tend
        v=state.reshape(6,n);der=rhs(t0,state).reshape(6,n)
        steady=float(np.max(abs(der)/np.maximum(abs(v),1)))
        if steady<2e-6:break
    c=speciate(v,cfg);f,psi,weights,se=fluxes(c,cfg,e,z,w,d)
    if steady>=2e-6:raise RuntimeError(f'Station did not reach steady balance: {steady}')
    if np.min(v[[0,1,3,4,5]]) < -1e-5:raise RuntimeError('Negative component in prescribed-current solve')
    if c.min() < -1e-5:raise RuntimeError('Negative species concentration, including Cs')
    # Balance includes axial convection, boundary transport and source terms.
    netcarbon=(f[0,0]-sum(f[:3,-1])-r['carbon']-np.dot(u/dx,((v[0]+v[1])-(upstream[0]+upstream[1]))*w))
    alkflux=W[2]@f
    netalk=(alkflux[0]-alkflux[-1]+r['electron']-np.dot(u/dx,(v[2]-upstream[2])*w))
    diag={'steady_scaled':steady,'carbon_balance_mol_m2_s':float(netcarbon),
      'alkalinity_balance_mol_m2_s':float(netalk),'max_charge_mM':float(np.max(abs(Z@c))),
      'max_water_relative_error':float(np.max(abs(c[3]*c[4]/k['KW']-1))),
      'max_ionic_drift_step':float(np.max(abs(Z[:,None]*psi[None,1:]))),
      'min_component_mM':float(np.min(v[[0,1,3,4,5]])),
      'gas_C_flux_mol_m2_s':float(f[0,0]),'liquid_edge_DIC_flux_mol_m2_s':float(sum(f[:3,cfg.n_catalyst])),
      'max_CO2_mM':float(c[0].max())}
    return {'components':v,'c':c,'flux':f,'psi':psi,'weights':weights,'electron_source':se,
      'z':z,'w':w,'d':d,'edges':e,'u':u,'dx':dx,'diagnostics':diag,'upstream':upstream}

def solve_channel(cfg,verbose=True):
    if not 0<cfg.chi<=cfg.porosity<=1:
        raise ValueError('Geometric-area Deff/Daq must be positive and no greater than liquid porosity in this model.')
    e,z,w,d=mesh(cfg);n=len(z);hyd=hydrodynamics(cfg);length=hyd['length_m']
    # Slight inlet refinement; full area-weighting uses actual dx.
    xe=length*np.linspace(0,1,cfg.n_axial+1)**1.3
    upstream=np.repeat(bulk(cfg)[:,None],n,axis=1)
    stations=[]
    for ix,dx in enumerate(np.diff(xe)):
        s=solve_station(cfg,upstream,dx);s['x_m']=.5*(xe[ix]+xe[ix+1]);stations.append(s)
        upstream=s['components']
        if verbose:print(json.dumps({'chi':cfg.chi,'x_step':ix+1,'n':cfg.n_axial,**s['diagnostics']}),flush=True)
    return stations

def summarize(cfg,stations):
    length=np.sqrt(cfg.area_m2)
    c=sum(s['c']*s['dx']/length for s in stations)
    nc=cfg.n_catalyst;cl=c[:,:nc]
    return {'chi':cfg.chi,'flow_mL_min':cfg.flow_mL_min,'n_catalyst':nc,'n_external':cfg.n_external,'n_axial':cfg.n_axial,
      'outer_domain_um':cfg.outer_domain_m*1e6,'carbon_flux_mmol_m2_s':1000*rates(cfg)['carbon'],
      'CL_CO2_min_mM':float(cl[0].min()),'CL_DIC_max_mM':float(cl[:3].sum(axis=0).max()),
      'CL_OH_mean_mM':float(cl[4].mean()),'CL_HCO3_mean_mM':float(cl[1].mean()),'CL_CO3_mean_mM':float(cl[2].mean()),
      'CL_pOH_min':float(-np.log10(cl[4].max()/1000)),'CL_pOH_max':float(-np.log10(cl[4].min()/1000)),
      'max_CO2_mM':float(c[0].max()),'max_ionic_strength_M':float(.5*(Z[:,None]**2*c).sum(axis=0).max()/1000),
      'max_steady_scaled':max(s['diagnostics']['steady_scaled'] for s in stations),
      'max_carbon_balance_mol_m2_s':max(abs(s['diagnostics']['carbon_balance_mol_m2_s']) for s in stations),
      'max_alkalinity_balance_mol_m2_s':max(abs(s['diagnostics']['alkalinity_balance_mol_m2_s']) for s in stations),
      'max_ionic_drift_step':max(s['diagnostics']['max_ionic_drift_step'] for s in stations)}

def isotope_channel(cfg,stations,gas_fraction,liquid_fraction):
    """Linear isotope transport through CO2 and the fast-proton-transfer DIC pool.

    Chemical carbon exchange remains finite in the catalyst AND the flowing
    liquid. A separate gas=1/liquid=0 solution identifies source contributions.
    Isotope reaction rates/diffusivities are equal here; this is no-KIE assumption.
    """
    k=constants(cfg);kh=cfg.hydration_s;ko=cfg.hydroxylation_M_s/1000
    upstream=np.repeat((bulk(cfg)[:2]*liquid_fraction)[:,None],len(stations[0]['z']),axis=1)
    answers=[]
    for s in stations:
        total=s['c'];n=len(s['z']);nc=cfg.n_catalyst;w=s['w'];d=s['d'][:3];psi=s['psi']
        btotal=total[1]+total[2]
        alpha=np.divide(total[1],btotal,out=np.ones(n),where=btotal>1e-30)
        eps=np.where(s['z']<cfg.catalyst_length_m,cfg.porosity,1.)
        af=eps*(kh+ko*total[4]);ar=eps*(kh/k['K1']*total[3]+ko*k['KW']/k['K1'])*alpha
        lam=rates(cfg)['carbon']*s['weights']/np.maximum(total[0],1e-30)
        cb=speciate(bulk(cfg)[:,None],cfg)[:3,0]*liquid_fraction
        g=1/(.5*w[:-1]/d[:,:-1]+.5*w[1:]/d[:,1:]);g=np.c_[g,d[:,-1]/(.5*w[-1])]
        def rhs(v,include_source):
            t=v.reshape(2,n);c=np.vstack([t[0],alpha*t[1],(1-alpha)*t[1]])
            f=np.zeros((3,n+1));f[0,0]=d[0,0]/(.5*w[0])*((cfg.CO2_star_mM*gas_fraction if include_source else 0)-t[0,0])
            cr=np.c_[c[:,1:],cb if include_source else np.zeros(3)]
            f[:,1:]=-g*(cr-c)-g*Z[:3,None]*.5*(cr+c)*psi[None,1:]
            dv=-(np.vstack([f[0,1:]-f[0,:-1],np.sum(f[1:,1:]-f[1:,:-1],axis=0)]))/w
            dv-=s['u']/s['dx']*(t-(upstream if include_source else 0))
            net=af*t[0]-ar*t[1]
            dv[0]-=net+lam*t[0];dv[1]+=net
            return dv.ravel()
        source=rhs(np.zeros(2*n),True);mat=lil_matrix((2*n,2*n))
        for j in range(2*n):
            b=np.zeros(2*n);b[j]=1;col=rhs(b,False)
            ii=np.flatnonzero(col);mat[ii,j]=col[ii,None]
        solution=spsolve(mat.tocsr(),-source).reshape(2,n)
        if solution.min() < -1e-6:raise RuntimeError('Negative isotope concentration')
        frac=np.divide(solution[0],total[0],out=np.zeros(n),where=total[0]>1e-20)
        product=float(np.dot(lam*solution[0],w)/max(rates(cfg)['carbon'],1e-30))
        answers.append({'components':solution,'fraction_CO2':frac,'product_fraction':product})
        upstream=solution
    return answers

def save(cfg,stations,name):
    OUT.mkdir(exist_ok=True,parents=True)
    length=np.sqrt(cfg.area_m2);c=sum(s['c']*s['dx']/length for s in stations);s=stations[-1]
    df=pd.DataFrame({'z_nm':s['z']*1e9,'liquid_y_um':(s['z']-cfg.catalyst_length_m)*1e6,
      'region':np.where(s['z']<cfg.catalyst_length_m,'catalyst','liquid')})
    for i,n in enumerate(SP):df[n+'_mM']=c[i]
    df['DIC_mM']=sum(c[:3]);df['pOH']=-np.log10(c[4]/1000);df['pH']=-np.log10(c[3]/1000)
    gas=isotope_channel(cfg,stations,1.,0.);liq=isotope_channel(cfg,stations,0.,1.)
    tg=sum(t['components']*s['dx']/length for t,s in zip(gas,stations))
    tl=sum(t['components']*s['dx']/length for t,s in zip(liq,stations))
    df['gas_fraction_CO2']=tg[0]/np.maximum(c[0],1e-30)
    df['gas_fraction_DIC']=sum(tg)/np.maximum(sum(c[:3]),1e-30)
    df['f13_CO2_13gas']=.0107+.9793*df['gas_fraction_CO2']
    df['f13_CO2_13liquid']=.0107+.9793*tl[0]/np.maximum(c[0],1e-30)
    df.to_csv(OUT/(name+'_mean_profile.csv'),index=False)
    frames=[]
    for ss,gg,ll in zip(stations,gas,liq):
        f=df[['z_nm','liquid_y_um','region']].copy();f['channel_x_mm']=ss['x_m']*1000;f['dx_m']=ss['dx']
        for i,n in enumerate(SP):f[n+'_mM']=ss['c'][i]
        f['gas_fraction_CO2']=gg['fraction_CO2'];f['liquid_fraction_CO2']=ll['fraction_CO2']
        frames.append(f)
    pd.concat(frames,ignore_index=True).to_csv(OUT/(name+'_axial_profiles.csv'),index=False)
    payload={'config':asdict(cfg),'hydrodynamics':hydrodynamics(cfg),'FE_rates':rates(cfg),'summary':summarize(cfg,stations),
      'diagnostics':[ss['diagnostics'] for ss in stations],
      'electrochemical_rate_fitted':False,'absolute_concentrations_validated':False,
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    payload['summary']['gas_source_fraction']=sum(t['product_fraction']*s['dx']/length for t,s in zip(gas,stations))
    payload['summary']['isotope_source_closure_abs_mM']=float(np.max(abs(tg+tl-np.vstack([c[0],c[1]+c[2]]))))
    payload['summary']['CL_gas_source_min']=float(df['gas_fraction_CO2'][:cfg.n_catalyst].min())
    payload['summary']['CL_gas_source_max']=float(df['gas_fraction_CO2'][:cfg.n_catalyst].max())
    (OUT/(name+'_manifest.json')).write_text(json.dumps(payload,indent=2),encoding='utf-8')
    return payload

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--chi',type=float,default=.25);p.add_argument('--nx',type=int,default=10)
    p.add_argument('--nc',type=int,default=24);p.add_argument('--ny',type=int,default=56)
    p.add_argument('--flow',type=float,default=10);p.add_argument('--outer-um',type=float,default=500)
    p.add_argument('--porosity',type=float,default=.5)
    p.add_argument('--co2-star',type=float,default=33.)
    p.add_argument('--pka1',type=float,default=6.35)
    p.add_argument('--allocation',choices=['CO2_weighted','uniform'],default='CO2_weighted')
    p.add_argument('--unresolved-fe',choices=['noncarbon','CO','formate'],default='noncarbon')
    p.add_argument('--name',default='current_central');a=p.parse_args()
    cfg=Config(chi=a.chi,porosity=a.porosity,n_axial=a.nx,n_catalyst=a.nc,n_external=a.ny,
      flow_mL_min=a.flow,outer_domain_m=a.outer_um*1e-6,CO2_star_mM=a.co2_star,pKa1=a.pka1,
      allocation=a.allocation,unresolved_FE_mode=a.unresolved_fe)
    result=solve_channel(cfg);print(json.dumps(save(cfg,result,a.name)['summary'],indent=2))
