# CO₂ reduction: spatial transport and isotope-source analysis

How can carbon supplied through a gas channel and carbon dissolved in the
electrolyte contribute differently to products on a porous Cu electrode?
This repository provides two complementary ways to examine that question:

1. **Analyze measured isotope signals.** Reprocess numerical NMR free-induction
   decays, fit isotope reporter patterns, analyze ethylene mass spectra, and
   convert product isotope fractions to gas-/liquid-carbon source fractions.
2. **Calculate spatial environments.** Solve CO₂, dissolved inorganic carbon
   (DIC), OH⁻, electrolyte-ion and isotope profiles across an idealized 200 nm
   aqueous-connected Cu path and its adjacent flowing liquid.

Experimental source fractions are observations to compare with transport
scenarios, not fitted transport parameters or unique reaction depths. All
coordinates, inputs and assumptions are defined here without paper figure numbers.

## Installation and quick start

Python 3.10 is the tested interpreter. From the repository root:

```sh
python -m pip install -r requirements.txt

# Experimental isotope analyses:
python -m isotope_analysis.nmr_templates --out runs/nmr
python -m isotope_analysis.gcms_analysis --out runs/gcms
python -m isotope_analysis.source_fractions --gcms runs/gcms/ethylene_response_envelopes.csv --out runs/sources

# Optional historical acetate-carbonyl sensitivity (see numerical limitations):
python -m isotope_analysis.acetate_carbonyl --out runs/acetate_carbonyl

# Saved spatial calculations and descriptive plots:
python validate_release.py
python plot_results.py --results results --out runs/spatial_profiles
python postprocess.py --results results --out runs/spatial_source_budgets
python profile_plots/code/plot_profiles.py
```

Use new output directories; archived inputs and reference results are not
overwritten. Plotting saved profiles does not rerun or refit the model.
PNG/SVG are default exports; PDF is optional.

## Experimental design and analytical outputs

The feeds are `12gas_12liquid` (natural-abundance reference), `13gas_12liquid`
(¹³CO₂ gas + ¹²CO₂ dissolved in liquid), and `12gas_13liquid` (¹²CO₂ gas +
¹³CO₂ dissolved in liquid). “¹²CO₂” retains natural ¹³C rather than meaning pure
carbon-12. The isotope experiments used a phosphine-modified Cu electrode at
300 mA cm⁻² in CO₂-saturated, acidified 0.5 M Cs₂SO₄ electrolyte, pH about 3.4.

| Analysis | Numerical inputs | Principal outputs |
|---|---|---|
| Proton NMR | Three complex FIDs and selected acquisition parameters | Spectra, reporter fits, component curves and residuals |
| Acetate carbonyl, exploratory | Liquid-labelled proton FID | 144 historical line-shape/window/coupling/baseline alternatives |
| Ethylene GC–MS | Four mass spectra, m/z 26–30 | Accepted response fits, isotope ranges and random-pairing comparisons |
| Carbon-source conversion | Declared NMR ranges and calculated ethylene ranges | Independent source estimates for each feed and their combined span |

See [analytical methods](docs/ISOTOPE_ANALYSIS.md) for equations, carbon positions
and uncertainty. Numerical FIDs and processed direct-¹³C spectra are included
without vendor operator/workstation metadata. The source dataset is
[archived on Figshare](https://doi.org/10.6084/m9.figshare.30215389).

Under 1.07% natural/99% enriched feed assumptions:

| Product (numerical reporter) | Gas-source fraction, span of both feed estimates |
|---|---:|
| Formate (carbon) | 93.9–99.1% |
| Acetate (methyl carbon) | 96.9–99.2% |
| Ethanol (methyl carbon) | 97.9–100% |
| n-Propanol (terminal methyl carbon) | 98.3–100% |
| Ethylene (mean carbon fraction) | 79.1–92.3% |

These are analytical-response/assignment spans, **not statistical confidence
intervals**. NMR endpoint ranges are curated analytical inputs; the preferred
fit does not automatically derive every endpoint. Raw and physically bounded
source conversions are both retained.

The optional acetate-carbonyl sensitivity has platform-dependent nonlinear
minima and linearized uncertainties. Its archived numerical output is retained
for inspection, but it is **not an input to the source-fraction table above**.
Structural checks do not establish convergence of those individual weak-signal
fits. See [validation scope](docs/VALIDATION.md#experimental-analysis-validation).

The separate [carbon-position ledger](isotope_analysis/inputs/carbon_position_assignments.csv)
also records predominantly gas-derived ethanol CH₂OH, both n-propanol methylenes,
and acetate carbonyl. Those assignments are retained without inventing precise
percentages from unreliable line-shape fits. CO remains qualitative because
nitrogen interference is unresolved; no quantitative CO fraction or bound is
reported. Methane signal was insufficient for isotope analysis.

## Spatial calculation: inputs and predictions

![Spatial profiles for the 10 µm CO₂ transfer-equivalent scenario](figures/spatial_profiles_d10.png)

* Prescribed: 300 mA cm⁻², measured mean product Faradaic efficiencies, a 200 nm
  liquid-connected path, reservoir chemistry and transport scenarios.
* Screened: CO₂ passive transfer-equivalent length δ = **10, 30, 50 µm** and
  effective diffusivity ratio χ = **0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 0.005,
  0.002, 0.001, 0.0005**.
* Solved: diffusion/electromigration, axial convection, finite reversible
  CO₂–DIC conversion, charge/speciation closure, current/carbon balance, and
  gas-/liquid-carbon transport through molecular CO₂ and ionic DIC.
* Assumed: the local carbon sink is proportional to molecular CO₂ and normalized
  to measured total carbon consumption. All carbon products share this spatial
  weighting. There is **no electrochemical rate calibration** against a chosen
  diffusivity or isotope result.
* Not independently predicted: product selectivities, molecular pathways,
  meniscus movement, pore wetting morphology or unique product-generation depths.

**δ is D/k, not a physical film boundary.** Species have different passive
transfer-equivalent lengths in one common 500 µm numerical liquid domain.
Enhanced effective shear defines scenarios, not actual pump velocity.
Smooth-wall geometry at 10 mL min⁻¹ gives about 108 µm; 10/30/50 µm are assumed
enhanced-transfer sensitivities, not measured values.

Read [model equations](docs/MODEL.md), [parameter sources](docs/PARAMETERS.md)
and [numerical validation](docs/VALIDATION.md) before interpreting concentrations.
Ideal activities and dilute-aqueous diffusivities are approximations in this
concentrated electrolyte. Saved cases passed recorded numerical gates; grid/domain
checks cover χ=0.05 only. Depleted-liquid qualifications and low-CO₂ data remain.

### Run new spatial cases

```sh
# Coarse smoke test, not a production-resolution prediction:
python model/species_transfer_model.py --delta 10 --chi 0.25 --nc 6 --ny 48 --nx 4 --out runs/smoke

# One production-resolution case:
python model/species_transfer_model.py --delta 10 --chi 0.25 --out runs/d10_chi0p25

# Thirty-case screen, optionally plus six verification cases:
python run_screen.py --out runs/new_screen --workers 2 --verification
```

The screen can take hours. Use `species_transfer_model.py` as the solver entry
point; its adjacent base modules are implementation dependencies.

## Files and reading guide

| Location | Contents |
|---|---|
| `isotope_analysis/` | Analytical code, numeric inputs and reference outputs |
| `docs/ISOTOPE_ANALYSIS.md` | Analytical equations, assumptions and input/output definitions |
| `model/` | Spatial solver and numerical contract checks |
| `fe_audit/faradaic_inputs.csv` | FE replicates, means and stoichiometric model inputs |
| `results/cases/` | 30 spatial scenarios and six verification cases, including raw numerical states |
| `profile_plots/` | Descriptive source tables, geometry illustration and editable plots |
| `figures/` | Preview spatial plots |
| `provenance/` | Parameter ledger, import hashes and release checksums |

The related [active-learning repository](https://github.com/Wang-Group/CO2RR-Phosphine-Modified-Cu)
contains catalyst-screening materials. See [release notes](docs/RELEASE_NOTES.md)
for version scope.

## Verification

```sh
python validate_release.py
python model/test_species_transfer.py
python -m unittest test_release_tools isotope_analysis.test_analysis isotope_analysis.test_nmr isotope_analysis.test_validation
```

Numerical reproducibility checks do not experimentally validate a particular
spatial mechanism.
