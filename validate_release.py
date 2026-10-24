"""Read-only verification of archived raw states, figure data and checksums."""
from pathlib import Path
from dataclasses import fields
import hashlib
import json
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'model'))
import species_transfer_model as s

def main():
    index=json.loads((ROOT/'results/case_index.json').read_text())
    assert len(index)==36 and sum(r['kind']=='screen' for r in index)==30
    profiles=pd.read_csv(ROOT/'profile_plots/data/profiles_all_transfer_scenarios.csv')
    count=0;maxima={}
    for row in index:
        folder=ROOT/'results/cases'/row['id']
        cfg=s.Config(**json.loads((folder/'config.json').read_text()))
        archived=json.loads((folder/'validated.json').read_text())
        assert archived['all_raw_restart_checks_passed']
        assert archived['raw_CO2_min_mM']>=0
        assert abs(archived['FE_rates']['carbon']-s.m.rates(cfg)['carbon'])<1e-14
        with np.load(folder/'raw.npz',allow_pickle=False) as arrays:
            for ix,dx in enumerate(arrays['dx_m']):
                state=arrays['components'][ix]
                assert state[0].min()>=0
                ok,recomputed=s.n.station_problem(cfg,arrays['upstream'][ix],float(dx))[1](state)
                assert ok,(row['id'],ix,recomputed['diagnostics'])
                for key in ['steady_scaled','carbon_balance_mol_m2_s','alkalinity_balance_mol_m2_s','max_charge_mM']:
                    maxima[key]=max(maxima.get(key,0),abs(recomputed['diagnostics'][key]))
                count+=1
        df=pd.read_csv(folder/'case_mean_profile.csv')
        assert np.isfinite(df.select_dtypes('number')).all().all()
        for species in s.m.SP:assert df[species+'_mM'].min()>=-1e-8
        assert len(df)==cfg.n_catalyst+cfg.n_external
        if row['kind']=='screen':
            f=profiles[np.isclose(profiles.delta_CO2_um,cfg.delta_CO2_um)&np.isclose(profiles.chi,cfg.chi)].reset_index(drop=True)
            assert len(f)==len(df)
            for column in ['z_nm','CO2_mM','DIC_mM','OH_mM','Cs_mM','f13_CO2_13gas','f13_CO2_13liquid']:
                assert np.allclose(df[column],f[column],rtol=2e-14,atol=1e-14),(row['id'],column)
        currents=pd.read_csv(folder/'current_and_tracer_axial.csv')
        cu=currents[currents.region.eq('catalyst')]
        integrated=np.sum(cu.total_current_mA_cm2_per_nm*(cfg.catalyst_length_m*1e9/cfg.n_catalyst)*cu.dx_m/np.sqrt(cfg.area_m2))
        assert np.isclose(integrated,cfg.current_mA_cm2,rtol=1e-9)
    manifest=ROOT/'provenance/release_checksums.json'
    checked=0
    if manifest.exists():
        for relative,expected in json.loads(manifest.read_text()).items():
            actual=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
            assert actual==expected,relative
            checked+=1
    print(json.dumps({'archived_cases_checked':len(index),'raw_station_states_rechecked':count,
      'production_figure_CSV_mapping':'passed','integrated_current':'passed',
      'raw_state_balance_gates':'passed','checksum_files_verified':checked,
      'maximum_recomputed_diagnostics':maxima,
      'scope':'Saved-state equation checks, not a full nonlinear campaign rerun or experimental validation'},indent=2))

if __name__=='__main__':main()
