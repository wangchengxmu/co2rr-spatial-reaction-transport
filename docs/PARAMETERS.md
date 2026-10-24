# Input provenance and numerical values

The machine-readable [parameter registry](../provenance/parameter_registry.csv)
distinguishes experimental, literature, derived, assumed and numerical inputs.
It retains the inherited scientific source audit; this release does not claim a
new verification of every primary numerical table. It has been updated to the
**actual corrected solver settings**, rather than older base-module defaults.

## Experimental constraints and explicit assumptions

* Current 300 mA cm⁻²; FEs are the means of three ordinary CO₂RR measurements,
  supplied in `fe_audit/faradaic_inputs.csv`. They are not complete FE balances
  measured separately in the isotope experiments.
* 200 nm path and 0.8 cm² exposed area map the experiment onto an ideal connected
  path and square planform. Liquid porosity 0.5 is assumed, not measured.
* Cell gap 2 mm and flow 10 mL min⁻¹ define the smooth-wall reference. Effective
  δ_CO₂=10/30/50 µm is a separate assumed sensitivity, not inferred from flow.
* Molecular CO₂=33 mM is a dilute-water boundary proxy, not a concentration
  assay in the working salt. Salting-out is omitted.
* Actual archived input pH is **3.4**; a reported experimental pH of 3.41 is not
  the exact number used in these runs. Activities are treated as concentrations.
* Natural/enriched 13C fractions 0.0107/0.99 are assumed source endmembers;
  no bottle-purity certificate or inlet assay is supplied here.

## Carbon kinetics

Forward baselines are `k_h=0.037 s⁻¹` and `k_o=12100 M⁻¹ s⁻¹`.
Inherited primary-source locator: Wang et al., *J. Phys. Chem. A* 2010,
114, 1734–1740, Table 1 p.1737, 25 °C; experimental conditions p.1735.
Reported uncertainties are ±0.002 s⁻¹ and ±400 M⁻¹ s⁻¹.
[DOI: 10.1021/jp909019u](https://doi.org/10.1021/jp909019u).
These are dilute-solution baselines, not measured concentrated-Cs coefficients.

With pK1=6.35 and pKw=14, reverse coefficients are **derived**:
`k_h/K1=82832.682 M⁻¹ s⁻¹` and `k_o Kw/K1=2.7088526×10⁻⁴ s⁻¹`.
They are not independent reverse measurements. Wang's reported reverse
hydroxylation value 4.00×10⁻⁴ s⁻¹ uses its own equilibrium closure; the inherited
source audit notes that it is constrained by microscopic reversibility.
The released parameter combination deliberately uses the selected ideal pK1,
not Wang's entire equilibrium/kinetic set. No instantaneous CO₂–DIC equilibrium
is imposed. Fast proton-transfer closures use pK2=10.33, pKS=1.99,
pKformic=3.75, pKacetic=4.756 and pKw=14.

## Diffusivities

All values below are `10⁻⁹ m² s⁻¹`; multiply by common χ inside Cu.

| Species | Value | Provenance/status |
|---|---:|---|
| CO₂ | 2.04527 | Zeebe temperature fit, 298.15 K |
| HCO₃⁻ | 1.10059 | Same primary fit |
| CO₃²⁻ | 0.804197 | Same primary fit |
| H⁺ | 9.2983996 | Dilute conductivity-derived benchmark |
| OH⁻ | 5.3054587 | Dilute conductivity-derived benchmark |
| Cs⁺ | 2.07 | Compiled dilute-ion proxy; primary numeric table unverified |
| HSO₄⁻ | 1.33 | Compiled dilute-ion proxy; primary numeric table unverified |
| SO₄²⁻ | 1.07 | Compiled dilute-ion proxy; primary numeric table unverified |
| HCOOH | 1.5 | Assumed neutral-acid proxy |
| HCOO⁻ | 1.454 | Compiled infinite-dilution value |
| CH₃COOH | 1.21 | Assumed neutral-acid proxy |
| CH₃COO⁻ | 1.089 | Compiled infinite-dilution value |

Zeebe, *Geochim. Cosmochim. Acta* 2011, 75, 2483–2498, Table 2 p.2491:
[10.1016/j.gca.2011.02.010](https://doi.org/10.1016/j.gca.2011.02.010).
Light et al., *Electrochem. Solid-State Lett.* 2005, 8, E16–E19, Table III:
[10.1149/1.1836121](https://doi.org/10.1149/1.1836121).
None of these is a measurement of transport inside this porous Cu in concentrated
Cs₂SO₄. Equal χ for all species, ideal activities, no ion pairing, no precipitation
and constant water activity are explicit simplifications.

## Lévêque and equilibrium-source qualifications

The 0.807 passive transfer coefficient is the rounded rectangular-wall similarity
constant. Source locator in the inherited audit: Sobolík et al., *J. Appl.
Electrochem.* 1998, 28, 329–335, Section 2.1 p.330 Eqs.(1a–b),
[10.1023/A:1003224018571](https://doi.org/10.1023/A:1003224018571).
Its unrelated probe calibration is not imported into this model.

Conventional pK values are not fitted or measured in the working electrolyte.
In particular Covington et al. report a bisulfate pK near 1.947, not exactly the
retained 1.99 ([10.1039/F19736900094](https://doi.org/10.1039/F19736900094)).
Dunn and Stokes provide context for ionization-dependent acid diffusion, not a
verified source for the exact neutral proxy values used here
([10.1071/CH9650285](https://doi.org/10.1071/CH9650285)).

## Corrected numerical controls

Production grid: 48 Cu cells / 96 liquid cells / 24 axial stations, outer domain
500 µm. Corrected BDF: `rtol=1e-9`, `atol=2e-10 mol m⁻³`; direct scaled root
polishing is also used. Required scaled steady residual <2e-6 s⁻¹. Proton-root
residual tolerance is 2e-12 mol m⁻³. Saved acceptance additionally checks carbon,
alkalinity, charge, water closure, positivity, drift-step size and source integrals.
These are implementation controls, not physical rates or experimental calibration.
