# Isotope signals, carbon positions and source fractions

## Experimental comparison

One isotope of CO₂ is supplied through the gas channel and the other dissolved
in electrolyte. The opposite arrangement is measured separately. Product isotope
composition reports its carbon reservoir, not by itself an elementary pathway
or unique reaction depth. Feed codes are defined in the root README and used
consistently throughout all tables.

## Proton NMR: numerical FID to reporter fraction

`inputs/nmr/*_fid.npz` contains 16,384 complex, unprocessed time-domain points
per acquisition. `acquisition.json` retains only fields needed for digital-filter
removal, spectral width and referencing, plus scan counts. The natural control
used 22 scans; each labelled sample used 256. Absolute intensities across samples
are not comparable without acquisition/concentration corrections.

Acquisition fields retain their Bruker names: `SW_h` is spectral width in Hz,
`SFO1` and processed reference `SF` are MHz, `O1` is Hz, and `GRPDLY` is digital
group delay in points. `DECIM`, `DSPFVS` and `AQ_mod` specify digital/acquisition
settings; `NS` counts scans. Keeping the processed `SF` value preserves the
original ppm reference, rather than inferring it from `O1` alone.

The preferred processing applies Bruker digital-filter removal (nmrglue 0.11),
0.3 Hz exponential broadening, zero-filling to 32,768 points, FFT/reversal, ACME
automatic phasing, and an asymmetric least-squares baseline (λ=2×10⁷, p=0.001,
12 iterations). Empirical natural-control multiplet templates represent ordinary
H–¹²C and symmetric H–¹³C branches. Component areas are nonnegative; constant
and linear local baseline coefficients are free.

```text
F13 = A13 / (A12 + A13)
A13 includes BOTH 13C-coupled branches.
```

Formate and acetate methyl have empirical-template and flexible symmetric-
satellite alternatives. Ethanol/n-propanol methyls are fitted jointly over
0.62–1.32 ppm. The ethanol branch near 1.251 ppm is not discarded as an arbitrary
nuisance; separate nuisance peaks occur near 1.1753 and 1.2623 ppm.

`nmr_templates` writes processed spectra, fitted/component curves and a summary.
The generated `analysis_settings.json` defines the numbered component columns.
The frozen reference summary supplies approximately:

| Reporter | Natural control | Gas-labelled | Liquid-labelled |
|---|---:|---:|---:|
| Formate empirical | 0.00788 | 0.93013 | 0.02031 |
| Formate flexible | 0.02488 | — | 0.04157 |
| Acetate methyl empirical | near zero | 0.95942 | 0.01902 |
| Acetate methyl flexible | 0.02473 | — | 0.02828 |
| Ethanol methyl joint | 0.01631 | 0.99908 | 0.02439 |
| n-Propanol methyl joint | near zero | about 1 | near zero |

Alternative fits are not independent experimental replicates. A natural-control
fit need not return exactly 0.0107; a near-zero fit does not prove absence.

### Curated ranges and independent positional assignments

`inputs/nmr_reporter_intervals.csv` contains **curated analytical ranges**, combining
template, integration and overlap/assignment judgments. They are declared inputs
to source conversion: a preferred-fit command does not generate every endpoint.
For example, the liquid-labelled ethanol range 0.012–0.031 includes alternatives
beyond the preferred joint fit. Acetate's gas-labelled value 0.96 is approximate;
equal endpoints do not imply zero measurement uncertainty.

The carbon-position ledger separately retains predominantly gas-derived ethanol
CH₂OH, both n-propanol methylenes, and acetate carbonyl. Position assignment is
stronger than the precision of exploratory automated line-shape fractions. An
overlap-sensitive numerical fit does not negate resolved isotope-dependent
patterns. Numerical reporter percentages are not complete molecular isotopologue
distributions.

Two processed direct-¹³C spectra and a peak-assignment list support carbon-region
inspection. They are not quantitative isotope-fraction measurements: relaxation,
NOE, overlap and concentration sensitivity were not calibrated for that purpose.
A missing automated peak-list entry alone is not a detection limit.

### Acetate carbonyl sensitivity

`acetate_carbonyl` uses the liquid-labelled FID, 0.25 Hz broadening and 65,536
points to test the central methyl doublet of ¹²CH₃–¹³COO⁻. Its methyl protons
couple across two bonds to carbonyl carbon. The 144 alternatives combine:

* Window halfwidths 0.018, 0.025, 0.035, 0.050 ppm.
* Two-bond coupling 5.5, 5.914 or 6.8 Hz.
* Gaussian, Lorentzian, pseudo-Voigt or Voigt lines.
* Constant, linear or quadratic baselines.

The 5.914 Hz center value derives from the gas-labelled four-line pattern.
Each model provides a best fraction and a formal linearized upper value. In
the lower-residual half, the median best fraction is about 0.00193 and maximum
formal upper about 0.0242. Residual ranking is not independent model validation.
These are line-shape sensitivity results, **not a calibrated detection limit
or a general carbonyl-label fraction**. The complete grid is retained.

## Ethylene GC–MS: fragments are not molecules

`inputs/ethylene_mass_spectra.csv` contains m/z 26–30 for an unlabelled standard,
natural reaction control and both isotope feeds. Exported masses were rounded
to nominal mass, retaining the maximum intensity if entries shared a mass.
The two controls have no exported m30 entry: original zero substitution is
retained and flagged by `entry_present_in_export`, not treated as a measured
zero. Inputs are sampled mass-spectrum intensities, not integrated isotope
chromatographic areas.

For nonnegative amounts `a=(a00,a01,a11)` of ¹²C–¹²C, mixed ¹²C–¹³C, and
¹³C–¹³C ethylene, a common fragment pattern is shifted by 0/1/2 mass units:

```text
    [k26   0    0  ]
    [k27  k26   0  ]
R = [k28  k27  k26]       y ≈ R a
    [ 0   k28  k27 ]
    [ 0    0   k28 ]
```

Natural-control probabilities are `((1-pN)², 2pN(1-pN), pN²)`, `pN=0.0107`.
Nonnegative least squares obtains reference fragment ratios. An 81×81 response
grid extends the two reference ratios by ±10%, with lower bound 0.25.
Acceptance requires relative L2 residual ≤ best+**0.02**. This is an assumed
sensitivity rule, not a statistical confidence criterion. Pure singly/doubly
labelled ethylene standards were unavailable; equally shifted fragmentation
patterns remain a model assumption.

The accepted families contain 644 and 374 fits. Mean ¹³C atom-fraction ranges
are 0.867584–0.914998 for gas-labelled feeding and 0.203920–0.215457 for
liquid-labelled feeding. All accepted rows are saved, preserving correlations.
Do not combine marginal endpoints from different fits into one distribution.

For each individual fit:

```text
p = P11 + P01/2
random mixed probability in one uniform independent pool = 2p(1-p)
mixed-pair deficit = 2p(1-p) - P01
```

The factor two includes both isotope orders. The same-fit deficit is positive
throughout both response families. This challenges one uniformly mixed independent
precursor pool; it does not uniquely identify a coupling intermediate. Spatial
or temporal heterogeneity and incomplete mixing remain possible.

CO is qualitative: N₂ contributes at m28 and m29, and no transferable nitrogen
subtraction or independently integrated CO isotope peaks were established. No
quantitative CO fraction, source bound or location is calculated. Methane signal
was insufficient for reliable isotope analysis.

## Isotope fraction to carbon source

For gas-source fraction `g`, assume two endmembers:

```text
F13_product = g F13_gas + (1-g) F13_liquid
g = (F13_product - F13_liquid)/(F13_gas - F13_liquid)
```

Defaults are natural ¹³C=0.0107 and enriched ¹³C=0.99, with equal isotope-dependent
rates. These are feed assumptions, not measured enrichments at each reacting
meniscus. The two feed estimates stay separate. Their combined span is the
smallest interval enclosing both estimates, not a statistical average assuming
identical local reservoirs.

`source_fractions` retains the raw algebraic value and the physically bounded
[0,1] presentation. Values above the assumed 99% isotope endmember yield raw
`g>1`; clipping is a declared bound, not measured perfect gas exclusivity.
Use `--natural` and `--enriched` for endmember sensitivity.

No unique spatial position follows from this arithmetic. Relating sources to
reaction regions additionally requires transport assumptions and consideration
of intermediate movement and product-specific activity.

## Provenance and reproduction boundaries

Measurements derive from the [source archive](https://doi.org/10.6084/m9.figshare.30215389).
`isotope_analysis/input_provenance.json` records numerical export transformations
and hashes. Only numerical arrays, allowlisted acquisition values and analytical
tables are deposited, not vendor operator/host/audit metadata or the entire
archive. Preferred fits, the acetate grid, ethylene families and source conversion
have distinct outputs to separate data, processing, fitted results, curated inputs
and interpretation.

Preferred-reporter regression uses a relative tolerance of 2×10⁻⁶ plus an
absolute floor of 2×10⁻¹⁰. Fitted positions and widths instead use a 10⁻⁵ ppm
absolute tolerance: 0.005 Hz at 500 MHz, about 1/50 of the zero-filled FFT
spacing. Different numerical libraries can return slightly different minimizers,
especially weak-control satellite positions. These computational tolerances are
not analytical confidence intervals. Frozen input and output file hashes still
require exact byte identity.
