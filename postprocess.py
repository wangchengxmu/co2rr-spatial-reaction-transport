"""Current-weighted source budgets and declared local-pairing null predictions.

No experimental isotope spectrum, response fitting, CO transport or kinetic fit
is performed. Source fields are read from archived spatial calculations.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=Path('results'))
    p.add_argument('--out',type=Path,default=Path('derived'))
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    summaries=[];pairs=[]
    index=json.loads((args.results/'case_index.json').read_text())
    for row in index:
        if row['kind']!='screen':continue
        folder=args.results/'cases'/row['id']
        cfg=json.loads((folder/'config.json').read_text())
        report=json.loads((folder/'validated.json').read_text())
        raw=pd.read_csv(folder/'case_axial_profiles.csv').query("region=='catalyst'").reset_index(drop=True)
        cur=pd.read_csv(folder/'current_and_tracer_axial.csv').query("region=='catalyst'").reset_index(drop=True)
        assert np.allclose(raw.channel_x_mm,cur.x_mm) and np.allclose(raw.z_nm,cur.z_nm)
        dz=cfg['catalyst_length_m']*1e9/cfg['n_catalyst'];L=np.sqrt(cfg['area_m2'])
        weight=cur.carbon_current_mA_cm2_per_nm.to_numpy()*dz*cur.dx_m.to_numpy()/L
        current=weight.sum();weight/=current
        expected=report['FE_rates']['carbon_electrons']*96485.33212/10
        assert np.isclose(current,expected,rtol=1e-9)
        gas=raw.gas_fraction_CO2.to_numpy();liquid=raw.liquid_fraction_CO2.to_numpy()
        assert np.max(abs(gas+liquid-1))<3e-6
        summaries.append(dict(case=row['id'],delta_CO2_um=cfg['delta_CO2_um'],
          reacting_carbon_gas_source_fraction=float(weight@gas),
          full_Cu_gas_source_min=float(gas.min()),full_Cu_gas_source_max=float(gas.max()),
          **report['summary']))
        for feed,f in [('13gas_12liquid',.0107+.9793*gas),('12gas_13liquid',.0107+.9793*liquid)]:
            mean=float(weight@f);variance=float(weight@((f-mean)**2))
            p00=float(weight@((1-f)**2));p01=float(weight@(2*f*(1-f)));p11=float(weight@(f*f))
            assert abs(p00+p01+p11-1)<1e-12
            assert abs(p01-(2*mean*(1-mean)-2*variance))<1e-12
            pairs.append(dict(case=row['id'],delta_CO2_um=cfg['delta_CO2_um'],chi=cfg['chi'],
              feed=feed,mean_13C_fraction=mean,P12_12=p00,P12_13=p01,P13_13=p11,
              single_fixed_pool_mixed=2*mean*(1-mean),variance_local_fraction=variance,
              assumption='local molecular-CO2 inheritance and common carbon-current weight'))
    assert len(summaries)==30 and len(pairs)==60
    pd.DataFrame(summaries).to_csv(args.out/'transport_summary.csv',index=False)
    pd.DataFrame(pairs).to_csv(args.out/'local_pairing_null_predictions.csv',index=False)
    print(json.dumps({'cases':len(summaries),'pair_predictions':len(pairs),'pair_identity_check':'passed',
      'no_parameter_fit':True,'no_low_CO2_masking':True,'not_product_specific_kinetics':True},indent=2))

if __name__=='__main__':main()
