# Model equations and assumptions

## Coordinates, variables and units

`z=0` is the gas-contact boundary of a connected pore-liquid path; Cu extends to
`Lc=200 nm`. External liquid begins immediately at `y=z−Lc=0`; there is no gap.
`x` follows flow over a square-window length `L=sqrt(0.8 cm²)=8.944 mm`.
The model does not assert that an entire real electrode is filled with liquid.
The gas-side GDL, pore geometry and meniscus movement are not resolved.

Concentrations are per liquid volume, in **mol m⁻³**, numerically equal to mM.
Fluxes are per geometric area. Catalyst liquid volume fraction ε is 0.5;
`D_eff,i=χ D_aq,i` there, and `D_eff,i=D_aq,i` outside. Homogeneous reaction
and storage terms are multiplied by ε in Cu. Thus χ is a geometric transport
factor, not the molecular diffusivity of liquid alone.

The 12 chemical species are CO₂, HCO₃⁻, CO₃²⁻, H⁺, OH⁻, Cs⁺, HSO₄⁻, SO₄²⁻,
HCOOH, HCOO⁻, CH₃COOH and CH₃COO⁻. Six components are transported:
molecular CO₂; ionic inorganic carbon B=HCO₃⁻+CO₃²⁻; total alkalinity A;
total sulfate S; total formate P; and total acetate Q.

`A = [HCO3] + 2[CO3] + [OH] − [H] + [SO4] + [HCOO] + [CH3COO]`.
Electroneutrality then implies `[Cs]=A+S`. Water, HCO₃⁻/CO₃²⁻,
HSO₄⁻/SO₄²⁻ and carboxylic-acid pairs are fast proton-transfer closures.
**Molecular CO₂ is not forced into local equilibrium with DIC.**

## Transport and external-transfer scenarios

For species i with charge qᵢ and dimensionless solution potential ψ=Fϕ/RT:

```text
N_i = −D_eff,i (∂c_i/∂z + q_i c_i ∂ψ/∂z)
ε ∂v/∂t = −∂(W N)/∂z − u ∂v/∂x + component sources
```

`W` is the explicit component matrix in `current_constrained_model.py`.
The solution-potential increment is solved at each face from the ionic current
needed to carry accumulated Faradaic charge. Electronic potential gradients,
double-layer structure and potential-dependent product kinetics are not solved.
The numerical station solve uses pseudo-time relaxation/direct-root polishing
to obtain steady balance, with an implicit upstream march in x. Axial diffusion
is neglected.

All species occupy the **same extended liquid domain**, normally 500 µm.
In external liquid `u(y)=γ_eff y`; in Cu `u=0`. The sensitivity closure is:

```text
γ_eff = D_CO2 L / (0.807 δ_CO2)^3
k_i = 0.807 (D_i² γ_eff / L)^(1/3)
δ_i = D_i/k_i = δ_CO2 (D_i/D_CO2)^(1/3)
```

These are passive Lévêque transfer-equivalent quantities. Reactive/migrative
profiles are solved explicitly and are not required to terminate at δᵢ.
The selected δ_CO₂ values are not fitted or measured. For 10 µm, the effective
shear is 34807.7 s⁻¹ versus the smooth-wall geometry reference 27.95 s⁻¹.
The legacy metadata field `mean_velocity_m_s` equals `γ_eff H/6`, an
**effective-shear-derived diagnostic, not the actual pump velocity**.
The actual mean Q/(width×gap) is 0.009317 m s⁻¹ at 10 mL min⁻¹.

## Finite reversible carbon conversion

In molar concentration convention the net CO₂-to-bicarbonate rate is:

```text
R = (k_h + k_o[OH−])[CO2]
    − ((k_h/K1)[H+] + k_o Kw/K1)[HCO3−]
```

`k_h=0.037 s⁻¹`, `k_o=12100 M⁻¹ s⁻¹`. Reverse coefficients follow from
the chosen ideal `pK1=6.35` and `pKw=14`. Bicarbonate/carbonate partition
uses `pK2=10.33`. Multiply R by ε in the geometric-volume component equations.
For mol m⁻³ units, second-order M⁻¹ s⁻¹ coefficients are divided by 1000.
This finite forward/reverse term operates **inside Cu and in external liquid**.

## Current, carbon consumption and alkalinity

For product p, measured FE is a fraction, nC,p carbon atoms and ne,p electrons:

```text
J_p = j FE_p/(ne,p F)
J_C = Σ nC,p J_p
w(z) = c_CO2(z) / ∫Cu c_CO2(z) dz   inside Cu; zero outside
r_C(z) = J_C w(z)
```

At **each axial station**, j=300 mA cm⁻²=3000 A m⁻².
The mean FE sum is 0.9568386333 without renormalization. The prescribed carbon
flux is 0.00533383047 mol C m⁻² s⁻¹. Carbon-product electron source follows w;
H₂ plus unassigned non-carbon current is uniform through Cu. Their sum integrates
to j/F=0.031092809 mol electrons m⁻² s⁻¹. Unassigned FE is 4.3161367%.

The total alkalinity source follows the complete electron source, not just
carbon reduction. Formate/acetate are added at their molecular rates and their
protonation-dependent anionic charge is already included in A. There is no
second subtraction of their charge. This source is not permanent free OH⁻:
reaction, buffering and transport partition alkalinity among the species.

No intrinsic electrochemical rate coefficient is fitted. Normalizing w enforces
measured total current/carbon consumption; it does not validate the assumed
spatial activity law. All carbon products use this same w, so differences in
measured product carbon sources cannot be independently predicted here.

## Boundary conditions and isotopes

At gas contact, molecular CO₂ is fixed to 33 mM; all other species have zero
normal flux. The channel inlet and outer-liquid boundary have CO₂=33 mM,
ideal concentration-index pH **3.4**, Cs=1 M and no organic products. Total
sulfate is derived from charge balance (about 0.5097245 M). The input liquid
is initially CO₂/DIC-equilibrated; local working profiles need not be.
There is no gas-side resistance, salting-out correction or finite-reservoir aging.

After the total concentration solution, solve two complementary linear source
problems through **CO₂ and ionic DIC**: gas=1/liquid=0 and gas=0/liquid=1.
They share the total-field migration potential, species partition and finite
forward/reverse carbon rates. Incoming bicarbonate can therefore diffuse,
regenerate CO₂, and supply a different location. The carbon sink for each source
uses its local molecular-CO₂ share. Isotope diffusivities/rates are assumed equal.

For gas-source fraction g and complementary liquid fraction l:

```text
13C fraction, 13CO2 gas + 12CO2 liquid = 0.0107 + 0.9793 g
13C fraction, 12CO2 gas + 13CO2 liquid = 0.0107 + 0.9793 l
```

Mean plotted percentages are **ratios of axially averaged labelled and total
concentrations**, not averages of percentages. Reaction-source budgets use
local carbon-consumption weighting before averaging. These two averages differ.
Low-CO₂ points are not masked. `postprocess.py` optionally evaluates local
random-pairing null probabilities, averaging p² and 2p(1−p) **after** calculating
local values; it does not simulate adsorbed CO or C–C bond formation.
