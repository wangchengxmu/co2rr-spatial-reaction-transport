"""Independent contract checks before screening the new transport closure."""
import json
import numpy as np
from dataclasses import replace
import species_transfer_model as s
m=s.m

def main():
    cfg=s.Config(n_catalyst=6,n_external=10,n_axial=8)
    e,z,w,d=m.mesh(cfg);N=len(z);up=np.repeat(m.bulk(cfg)[:,None],N,axis=1)
    # A physically speciable, nonuniform state, with positive product pools.
    v=up.copy();v[0]*=np.linspace(.4,1,N);v[1]+=1
    v[2]+=15*np.exp(-z/30e-6);v[4]+=.05;v[5]+=.02
    rhs,assess=s.n.station_problem(cfg,up,.001)
    jac=s.jacobian(rhs,0,v.ravel(),cfg.n_catalyst).toarray()
    base=rhs(0,v.ravel());full=np.zeros_like(jac)
    for col in range(len(base)):
        vv=v.ravel().copy();step=1e-6*max(abs(vv[col]),1);vv[col]+=step
        full[:,col]=(rhs(0,vv)-base)/step
    error=float(np.max(abs(jac-full)/np.maximum(abs(full),1)))
    assert error<1e-9,error
    c=m.speciate(v,cfg);weights,se=m.spatial_sources(c,cfg,w);rates=m.rates(cfg)
    assert abs(weights@w-1)<1e-12
    assert abs(se@w-rates['electron'])<1e-12
    assert np.ptp(weights[:cfg.n_catalyst])>0
    assert np.allclose(se[:cfg.n_catalyst]-rates['carbon_electrons']*weights[:cfg.n_catalyst],rates['noncarbon_electrons']/cfg.catalyst_length_m)
    passive=[]
    for delta in [10,30,50]:
        c=replace(s.Config(),delta_CO2_um=delta)
        hy=s.hydro(c)
        lengths=np.array(list(hy['delta_species_m'].values()))
        assert np.allclose(lengths/(delta*1e-6),(m.D/m.D[0])**(1/3))
        for nx,ny in [(24,96),(96,192)]:
            row=s.passive(c,nx=nx,ny=ny);passive.append(row)
            assert abs(row['relative_k_error'])<.06,row
        # The D^(2/3) scaling is checked independently for four diffusivities.
        for idx in [1,2,3,4]:
            row=s.passive(c,diffusivity=m.D[idx],nx=96,ny=192);passive.append(row)
            assert abs(row['relative_k_error'])<.025,row
    result={'passed':True,'jacobian_max_scaled_difference':error,'source_integrals_passed':True,
      'passive_benchmarks':passive,'definition':'Numerical/implementation checks, not validation against experimental concentration or isotope data.'}
    s.write(s.ROOT/'contract_tests.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
