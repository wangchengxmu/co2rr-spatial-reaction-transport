# Verification status and limitations

## Archived results

The release contains the complete 30-case production screen and six verification
cases. All passed the recorded numerical acceptance gates. Original source-code
and FE-input SHA256 values matched the campaign manifest at release preparation.
Raw states (`raw.npz`, loaded with `allow_pickle=False`) and upstream fields are
retained, allowing the balance equations to be reevaluated rather than trusting
plot CSVs or a solver-success flag alone. Concentrations are not rounded or
clipped during packaging.

The six additional cases comprise doubled Cu/liquid/axial meshes and 1000 µm
outer domains at χ=0.05 for each δ_CO₂. These **do not establish mesh convergence
of every χ**, especially sharply depleted low-diffusivity cases. Detailed changes
are in `results/verification_comparisons.json`.

Contract checks independently compare the colored/nonlocal Jacobian with full
finite differences, integrate the carbon/current source allocation, and compare
passive finite-volume transfer with the analytical Lévêque scaling. Refined
passive CO₂ equivalent lengths are approximately 10.006, 30.018 and 50.029 µm.
This verifies implementation, not real-electrode morphology or enhanced transfer.

## Depleted-liquid isotope qualification

Low-CO₂ values are retained. At δ_CO₂=50 µm, χ=0.001 and 0.002, some depleted
external-liquid cells have relative molecular-CO₂ tracer-reconstruction errors
of about 1.59% and 1.30% at total CO₂ near 1.4×10⁻⁸ and 1.6×10⁻⁸ mM. They pass
the small absolute-error gate, but percentages in these cells need further
numerical verification. The complete region/case audit is supplied in
`results/Cu_isotope_distribution_audit.json` and `.csv`. No masking removes them.
The maximum relative source-reconstruction error reported within Cu is about
2.2×10⁻⁶. Numerical consistency alone does not prove physical accuracy.

## Interpretation boundaries

* Ideal solution chemistry at formal ionic strength around 1.5 M is an
  approximation. Large DIC or free-OH⁻ values are not validated concentration
  predictions. DIC is not limited by molecular CO₂ solubility, but must satisfy
  flux/charge/alkalinity balance; that is not a proof that omitted activities,
  association, salting-out or precipitation are negligible.
* Effective shear is a sensitivity device. It is not a measured pump speed or a
  simulation of moving menisci. The physical plausibility of 10 µm transfer
  cannot be established by the isotope result used to interpret that scenario.
* All carbon-product currents share a common CO₂-dependent spatial weight.
  Product-resolved isotope data are observations to compare against this model,
  not independent model-predicted selectivities or unique depth assignments.
* No mobile CO pool, sequential CO₂→CO→C₂ chemistry, surface intermediate,
  isotope-dependent kinetics, transient wetting or precipitation is simulated.
* Local pairing is optional post-processing with a declared precursor-inheritance
  assumption. A mixed-pair deficit can arise from spatial/temporal heterogeneity;
  it does not by itself identify a single elementary coupling pathway.

## Commands

`python validate_release.py` checks all archived raw states against the corrected
equations and input rates, the 30-case figure CSV mapping, current integrals,
positivity and a source/archived-data checksum manifest (regenerable preview
figures are excluded from the checksum manifest). It does not rerun every nonlinear
solve. `python model/test_species_transfer.py` runs the Jacobian/source/passive
benchmark tests. A separate smoke run demonstrates that the packaged solver can
produce a new result. Full-screen recalculation is an explicit separate command.

## Experimental analysis validation

Preferred proton-reporter fits are regression-tested with unit-aware numerical
tolerances. Ethylene response families and source conversion have independent
tests, including natural-abundance correction, both feed directions and the
factor two in random mixed-pair probability. Frozen numerical input/output
hashes remain byte-exact checks.

The optional historical acetate-carbonyl grid is different: strict regression
of all fitted nuisance parameters fails in the checked Linux environment.
Some individual minima and residuals differ materially, not just in their last
digits. The original algorithm and Windows snapshot are retained with this
limitation; they are not inputs to the accepted product source-fraction table.
CI uses explicit `--contract-only` validation for this grid: all 144 setting
combinations, finite values, declared bounds and summaries recomputed from the
generated table must agree internally. Strict reference differences are still
reported in full and retained with the numerical tables as a CI artifact.
A passing structural contract must not be described as convergence validation
or exact reproduction of these exploratory fits.
