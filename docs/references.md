# References, and what was actually read

Every source cited in this documentation, with a note on whether it was read
directly or quoted through someone else. That distinction matters: a value
carried through two papers can pick up a change of definition on the way, and
a citation implies the reader can check it against the same thing the author
saw.

## Read directly

**Spicer, T.O. and Havens, J.A. (1989).** *User's Guide for the DEGADIS 2.1
Dense Gas Dispersion Model.* EPA-450/4-89-019. The model this package ports;
the Fortran source and the five test cases are vendored in `reference/`.

**Zapert, J.G., Londergan, R.J. and Thistle, H. (1991).** *Evaluation of Dense
Gas Simulation Models.* EPA-450/4-90-018. Read for Table 5-7 (fractional bias
by model) and sections 4.2.1 and 4.3 (input and output assumptions). The
statement that measurements were taken at 1 m and compared against ground-level
predictions is section 4.3.

**AEA Technology (2001).** *Integral Modelling of the Dilution and Lift-off of
Ground Based Buoyant Plumes and Comparison with Wind Tunnel Data.*
AEAT/NOIL/27328006/001, EC URAHFREP Work Package 7. The lift-off model in
`addons/liftoff.py` is implemented from this.

**AEA Technology (2001).** *Model Predictions Compared with URAHFREP Campaign
2 Field Trial Data.* Read for the lift-off parameter definitions and for the
finding that peak or short-averaged data agrees best.

**Witcofski, R.D. and Chirivella, J.E. (1984).** "Experimental and analytical
analyses of the mechanisms governing the dispersion of flammable clouds formed
by liquid hydrogen spills." *Int. J. Hydrogen Energy* **9**(5), 425–435.
Tables 1 to 4 are used directly.

**Middha, P., Ichard, M. and Arntzen, B.J.** *Validation of CFD Modelling of
LH2 Spread and Evaporation Against Large-Scale Spill Experiments.* GexCon.
Read for the stability-class sensitivity and the Test 6 source conditions.

**PRESLHY D3.2** and **D3.6** (Fuel Cells and Hydrogen Joint Undertaking,
grant 779613). D3.6 supplies the E3.5 conditions, the sensor table and the
instrument ranges; D3.2 the integral-model comparisons and the 72 K validity
limit on the adiabatic mixing method.

**Lyons, K., Coldrick, S. and Atkinson, G. (2023).** *Fuel Cells and Hydrogen
Joint Undertaking (FCH JU); Summary of experiment series E3.5*, Karlsruhe
Institute of Technology, DOI 10.35097/1481.  The public archive catalogue was
read directly for its instrument inventory and archive description.  It lists
the H2/O2/thermal/pressure/mass-flow/weather/video channels and reports
condensed air near the release, but contains no listed particle-resolved or
separate gas/particle-velocity channel:
<https://doi.org/10.35097/1481>.

**Norwegian Defence Research Establishment (2021).** *Large scale leakage of
liquid hydrogen (LH2) — tests related to bunkering and maritime use.* FFI
Report 20/03101. Pages 15--17 were visually checked: seven outdoor tests use
pad/field thermocouples, oxygen-depletion gas measurements, and 5/10 m mast
weather readings; page 17 identifies the 30/50/100 m, 0.1/1.0/1.8 m field
array. No particle-size, phase-inventory or gas/particle-velocity instrument
is described. The local source PDF is not redistributed.

## Primary-source correction

**Papanicolaou, P.N. & List, E.J. (1988).** "Investigations of round vertical
turbulent buoyant jets." *Journal of Fluid Mechanics* **195**, 341–391.
doi:10.1017/S0022112088002447. Open-access record: CaltechAUTHORS
`k4bya-erm31`.

The primary paper has now been read. Section 3.1, pp. 353–354, derives a jet
entrainment coefficient of **0.0545**, a plume coefficient of **0.0875**, and
a pure-plume Richardson number of **0.716** from the experiment. The often
quoted 0.0533, 0.0833 and `Ri_p = 0.557` are Fischer et al. (1979)'s proposed
values, quoted by Papanicolaou and List for comparison; their earlier
attribution through two secondary sources is superseded.

This changes the evidentiary status in both directions. The coefficient is
derived from a primary measurement, not merely corroborated, but the adopted
number must be 0.0875 rather than 0.0833. `hydrogen_jet` now uses the measured
value. Historical parameter-dump tests request Fischer's value explicitly so
that reconstruction of the old run remains exact instead of silently freezing
the model at the old coefficient. HyRAM's 0.082 is independent corroboration,
not the derivation.

**List, E.J. (1982).** "Turbulent jets and plumes." *Annu. Rev. Fluid Mech.*
**14**, 189–212. doi:10.1146/annurev.fl.14.010182.001201. This remains a valid
background reference, but it is no longer the source used for the adopted
coefficient.

## Corrections adopted on a relation this work has not verified at source

The same test applies to the other two corrections, and one of them does not
pass cleanly either.

**Ricou, F.P. and Spalding, D.B. (1961).** "Measurements of entrainment by
axisymmetrical turbulent jets." *J. Fluid Mech.* **11**, 21–32.

This is the basis for scaling the shear entrainment as
`sqrt(rho_ambient / rho_jet)`. The citation is standard and the relation is in
routine use, including by Panda and Hecht for cryogenic hydrogen. But the
density-ratio scaling that everyone attributes to these experiments has
recently been challenged: an analysis in the same literature argues that the
1961 measurements do not by themselves settle the dependence of the
entrainment coefficient on the density ratio, and that the physical mechanism
is therefore not established by them.

This does not overturn the correction — it is applied here in the direction the
relation states. The former 0.70-to-0.76 width comparison used a measured
e-folding width under a sigma label; on the common standard-deviation
definition it is about 0.99-to-1.08. It does mean
the correction should be described as **a widely used scaling whose
experimental basis has been questioned**, rather than as a measured law. The
sign was inverted on first attempt here, which is itself a reminder that the
relation was being applied from memory of its form rather than from the paper.

**The expanded-source geometry correction is clean.** Its basis is
EPA-450/4-90-018 §4.3, which is in hand and quoted directly: source density,
composition and area must describe the same fully expanded plane. The later
air-entrainment step is a separate boundary. DEGALI now conserves total
mixture momentum there with the no-slip Sandia balance rather than extending
the earlier density-only velocity scaling across entrained mass.

**Papanicolaou and List (1988)** has now been read at source. **Kaminski, Tait
and Carazzo (2005)** remains a secondary comparison for this project rather
than the basis of the adopted value.

**Xiao, J. and co-workers (2009).** "Non-Boussinesq integral model for
horizontal turbulent buoyant round jets." *Sci. Technol. Nucl. Install.*
doi:10.1155/2009/862934. Abstract and quoted passages read; the full
derivation has not been.

**Mack, A. and co-workers (2023).** *Process Safety and Environmental
Protection* **176**. Read for the pressure-drag finding and the Witcofski
conditions table; the EFFECTS implementation itself has not been examined.

**Giannissi, S. and co-workers.** ADREA-HF simulations of cryogenic hydrogen
releases. Read for the finding that humidity condensation dominates LH₂ cloud
buoyancy and that air condensation is confined near the release.

**Li, X., Zhang, J., Wang, Y., Christopher, D.M., Yin, Q. and An, G. (2026).**
"Modeling of cryogenic hydrogen jets with air condensation."
*International Journal of Hydrogen Energy* **250**, article 156128.
doi:10.1016/j.ijhydene.2026.156128. The final article was supplied and read
directly. Equations 13--24 and Fig. 4 are reproduced in
`addons.cryogenic_air`; the phase-domain and energy-balance audit is in
`docs/li2026-air-condensation-audit.md`.

**Barrett, C.S., Meyer, L., Greer, S.C. and Wasserman, J. (1968); Erratum
(1972).** "Nitrogen--Oxygen Phase Diagram," *Journal of Chemical Physics*
**48**(6), 2670--2673, doi:10.1063/1.1669500; erratum **56**, 681,
doi:10.1063/1.1676933.  The author-hosted public copy was checked on
2026-09-15 and contains the erratum, not a usable original phase-diagram
table. It records a Gibbs-phase-rule correction at 43 K and explicitly says
the corrected region had no experimental data. Search metadata describes the
original X-ray observations only over 21--50 K. This source therefore
supports a low-temperature mixed-solid warning, but does not supply the
59--63 K mixture thermodynamic closure; see
`cold-mixed-phase-evidence.md`.

**Kochenburger, T.M., Grohmann, S. and Oelrich, L.R. (2015).** "Evaluation
of a two-stage mixed refrigerant cascade for HTS cooling below 60 K."
*Physics Procedia* **67**, 227--232, doi:10.1016/j.phpro.2015.06.039.
The open article's Fig. 1 renders experimental N2--O2 liquidus and solidus
points over roughly 48--64 K, explicitly attributing the phase diagram to
Barrett et al. and the experimental points to Ruhemann et al. It is an
independent public check that oxygen depresses nitrogen freezing and that the
59--63 K local states can lie in a mixed-phase boundary region. It does **not**
provide a tabulated, pressure-qualified mixture Gibbs/enthalpy closure, so it
is used only to state the applicability boundary and not to create numerical
DEGALI coefficients; see `cold-mixed-phase-evidence.md`.

**Ruhemann, M. (1936).** "Equilibrium Curves of Low-Melting Mixtures,"
*Physikalische Zeitschrift der Sowjetunion* **16**, 67--82 (special issue,
June). The accessible NBS low-temperature oxygen bibliography identifies this
as experimental two-phase-equilibrium and specific-heat work on O2--N2 over
50--80 K, with graphs and apparatus. Its numerical paper is not publicly
available in the sources checked, so the bibliographic description is not
converted into a Cp, enthalpy, phase-fraction, or mixture-EOS coefficient.
It is the precisely identified primary-source target for the 59--63 K
mixed-phase energy closure; see `cold-mixed-phase-evidence.md`.

**Zhang, J.X., Ba, Q.X., Xiao, J.S., Christopher, D.M., Liu, Y., Yao, C.Y. and
Li, X.F. (2023).** "Analytical model of cryogenic hydrogen releases." ICHS
conference paper 189. Read directly from the public HySafe PDF:
<https://hysafe.info/uploads/papers/2023/189.pdf>. Its Zone-II
model prints mass, momentum and energy balances for entrained air and condensed
nitrogen and is the public precursor to Li et al. (2026). It cannot safely be
implemented verbatim: the printed latent-heat term lacks a mass-flow factor,
one energy expression substitutes mixture density for enthalpy, and the
condensed-nitrogen interpolation has a reversed sign in the phase-change
interval. The final 2026 paper resolves the intended printed form, while the
direct audit finds separate phase-domain and energy-bookkeeping limitations.

**Sandia National Laboratories (2025).** *Hydrogen Plus Other Alternative
Fuels Risk Assessment Models (HyRAM+) Version 6.0 Technical Reference Manual*.
SAND2025-04942, section 3.2.3, equations 58--64.  The manual states the mass,
momentum and energy balances for the optional initial entrainment/heating
plug-flow zone and records that its 0 K default leaves it disabled in the GUI:
<https://www.osti.gov/servlets/purl/2563814>.  The official PDF is retained in
`reference/lh2/` with a checksum.

**Sandia National Laboratories (2023).** *Technical Justifications for Liquid
Hydrogen Exposure Distances*, SAND2023-12548.  Read directly, especially
section 3.  It documents the HEM/MLM source-flow comparison and its limited
DNV/PRESLHY concentration comparisons.  Crucially, it states that its
dispersion and flame models are gas-only after the orifice and neglect air or
humidity condensation.  It is therefore a useful independently published
scope boundary, not evidence that a gas-only calculation and DEGALI's optional
condensed-air branches have been compared under identical assumptions:
<https://energy.sandia.gov/download/80280/>.

**Sandia National Laboratories, HyRAM jet model.** The current `_dev_plug`
implementation was read directly at commit
`b45abf9a6d995951311be6aad836f1874e4d420b`. Its initial
entrainment/heating step uses `v_out = v_in * mdot_H2_in / mdot_out`, the
no-slip momentum balance adopted for DEGALI's corrected equivalent-source
boundary. It also confirms that the minimum handoff temperature remains a
user input. This supports the balance structure, not a pressure-thrust value
for the preceding under-expanded zone:
<https://github.com/sandialabs/hyram/blob/master/src/hyram/phys/_jet.py>.

**Hecht, E.S. and Panda, P.P. (2019).** "Mixing and warming of cryogenic
hydrogen releases." Controlled 1 and 1.25 mm cryogenic hydrogen jets at
2--5 bar absolute. The measured Gaussian mass-fraction half-width growth was
0.07069 in the final manuscript (0.06503 in the 2017 conference version),
below the roughly 0.10--0.11 room-temperature values cited by the
authors. This directly argues against repairing DEGALI by arbitrarily
increasing total entrainment. Read from Sandia's publication record and the
official OSTI manuscript: <https://www.sandia.gov/research/publications/details/mixing-and-warming-of-cryogenic-hydrogen-releases-2019-04-02/>,
<https://www.osti.gov/servlets/purl/1529288>. The aggregate figure legends
also contain an untabulated 4 bar/45 K series, so the nine-condition model
comparison is marked provisional pending fit provenance.

**Hecht, E.S. and co-workers (2023).** "Concentration fluctuations and
flammability of cryo-compressed hydrogen and methane jets." *Fuel* **358**,
130230. doi:10.1016/j.fuel.2023.130230. The OSTI accepted manuscript was read
directly. Its 1.25 mm, 3 bar, 51 K hydrogen case provides concentration RMS,
intermittency and radial PDF-shape evidence, not velocity TKE or epsilon:
<https://www.osti.gov/biblio/2311268>.

**Mohammadpour, J. and co-workers (2024).** "An efficient dimensionality
reduction approach for modelling cryogenic hydrogen release."
*International Journal of Hydrogen Energy* **91**, 649--658.
doi:10.1016/j.ijhydene.2024.10.182. The open institutional manuscript was
read directly. Its Figure 3 mean and axial-RMS velocity contours for a
5 bar, 50 K, 1 mm LES are used only for a positive-semidefinite TKE screening
envelope, because the paper explicitly notes that the compared Sandia report
did not publish velocity distributions:
<https://research-management.mq.edu.au/ws/portalfiles/portal/399409510/394825139.pdf>.

**Hecht, E.S. and Panda, P.P. (2018).** *Experimental validation of a model
for cryogenic hydrogen jet dispersion*, SAND2018-2834C. The public Sandia
presentation was read directly through OSTI record 1503433. It confirms
simultaneous Raman imaging and condensed-moisture-particle PIV, plus displayed
mean-velocity panels, but provides no downloadable velocity-RMS, Reynolds
stress, epsilon, integral-scale or tracer-response data. It is therefore
evidence of the measurement capability and a provenance caution, not a TKE
calibration source: <https://www.osti.gov/biblio/1503433>.

**Saini, D. and co-workers (2024).** "Effect of storage conditions on the
characteristics of cryogenic hydrogen jet dispersion." *International Journal
of Hydrogen Energy* **67**, 1--15. doi:10.1016/j.ijhydene.2024.04.040. This
open cryogenic 3--7 bar, 50--70 K RANS study documents a conventional
`k-epsilon` comparison formulation (`C_mu=.09`, `sigma_k=1`), but does not
provide a matched experimental epsilon or integral-length-scale profile. It
is therefore evidence for transparent notation only, not a DEGALI TKE closure:
<https://doi.org/10.1016/j.ijhydene.2024.04.040>.

**Schiller--Naumann finite-Re isolated-sphere drag relation.** Its
`24/Re_p * (1 + 0.15 Re_p**0.687)` form and the conventional `Re_p <= 1000`
range were checked against the open *Journal of Fluid Mechanics* discussion
of finite-Re point-particle drag:
<https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/new-paradigm-for-computing-hydrodynamic-forces-on-particles-in-eulerlagrange-pointparticle-simulations/C653BA1C01108A86F99CF0F921931B44>.
DEGALI uses only the local response-time form behind an explicit opt-in
helper; it does not use this incompressible isolated-sphere relation as an
LH2 particle-size or drag-law selection.

**Davies, C.N. (1945), Cunningham slip correction.** The conventional
`C_c = 1 + Kn[1.257 + 0.400 exp(-1.10/Kn)]` form and the diameter convention
`Kn=2 lambda/d` were checked against the U.S. EPA's reproduction of Davies'
equation and definition:
<https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=9101ZF9P.TXT>. It is exposed
only as a declared-input low-Re rarefaction screen, not as a cryogenic-mixture
mean-free-path or combined finite-Re drag closure.

**Mohammadpour, J., Ba, Q., Li, X. and Salehi, F. (2026).** "Innovative
approaches for predicting cryogenic hydrogen behaviour." *International
Journal of Heat and Fluid Flow* **117**, 110025.
doi:10.1016/j.ijheatfluidflow.2025.110025. This open article uses 540 LES
snapshots of the Sandia 5 bar/50 K jet for a POD--BiLSTM reconstruction, but
states that the numerical data are available on request. Its snapshot results
are not treated as public experimental RMS data or an LH2 calibration target:
<https://doi.org/10.1016/j.ijheatfluidflow.2025.110025>.

**Li, X., Zhang, J., Wang, Y., Christopher, D. M. and Yin, Q. (2026).**
"Modeling of cryogenic hydrogen jets with air condensation." *International
Journal of Hydrogen Energy* **250**, 156128.
doi:10.1016/j.ijhydene.2026.156128. Its initial entrainment/heating zone
adds equilibrium nitrogen-condensation latent heat and a stationary-condensate
momentum term. DEGALI independently implements this as an opt-in limiting
comparison, retaining its plug-flow, nitrogen-only and equilibrium limits:
<https://doi.org/10.1016/j.ijhydene.2026.156128>.

**Sun, R., Pu, L., He, Y., Wang, T. and Tan, H. (2024).** "Phase change
modeling of air at the liquid hydrogen release." *International Journal of
Hydrogen Energy* **50**, 717–731. doi:10.1016/j.ijhydene.2023.06.201. The
abstract and indexed summary were read. Condensed air increases mixture
density, while latent heat raises temperature and buoyancy; the reported
continuous turbulent release is only weakly sensitive to the air-condensation
coefficient.

**National Bureau of Standards (1955).** *Tables of Thermal Properties of
Gases*, Circular 564, Table 7-11.  The solid-N2 vapour-pressure table supports
the sub-triple-point nitrogen correlation used by `air_saturation_pressure`:
<https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular564.pdf>.

**Georgia Institute of Technology, Project A-593.** Table 35 transcribes
Aoyama and Kanda's measured solid-O2 vapour pressures from 36 to 54.36 K; the
monotone points support the oxygen branch of `air_saturation_pressure`:
<https://repository.gatech.edu/bitstreams/3b16b36e-0022-4035-b21f-a0bba98886b4/download>.

**NIST Chemistry WebBook, SRD 69: Argon.** Read for atmospheric-argon phase
screening: molecular weight, approximately 83.8 K triple temperature and
approximately 0.689 bar triple pressure. The implemented sub-triple branch is
a declared constant-latent Clausius--Clapeyron continuation, not a hidden
CoolProp extrapolation:
<https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440371&Mask=24>.

**Murphy, D.M. and Koop, T. (2005).** "Review of the vapour pressures of ice
and supercooled water for atmospheric applications." *Quarterly Journal of
the Royal Meteorological Society* **131**, 1539--1565.
doi:10.1256/qj.04.94. Equation 7 supplies the stable H2O ice--vapour pressure
used above 110 K by the optional humid equilibrium-phase table. Below 110 K
the code labels its boundary-anchored constant-latent continuation as an
extrapolation rather than experimental frost kinetics:
<https://doi.org/10.1256/qj.04.94>.

**IAPWS R10-06 (2009).** *Revised Release on the Equation of State 2006 for
Ice Ih.* Table 6 gives 916.7094922 kg/m3 at the water triple point; that
value is the constant solid-H2O volume bound used by the humid equilibrium
option, rather than liquid-water density:
<https://www.iapws.org/relguide/Ice-Rev2009.pdf>.

**Luna, A.J. and co-workers (2018).** "Measurements of Enthalpy of
Sublimation of Ne, N2, O2, Ar, CO2, Kr, Xe, and H2O Using a Double Paddle
Oscillator." Read from the official NIST copy for the 7.79 kJ/mol argon
sublimation value used only in the pre-registered sensitivity:
<https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=921657>.

**2022 CODATA adjustment.** The exact Stefan--Boltzmann constant
`5.670374419e-8 W m-2 K-4` supports the deliberately perfect-black radiation
upper bound:
<https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf>.

**Leachman, J.W., Jacobsen, R.T., Penoncello, S.G. and Lemmon, E.W. (2009).**
"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and
Orthohydrogen." *Journal of Physical and Chemical Reference Data* **38**,
721--748. doi:10.1063/1.3160306. The NIST equations and their CoolProp
implementations supply the separate fuel caloric tables used in the
pre-registered spin-isomer sensitivity. The para-hydrogen candidate is
numerically conservative but worsens the Raman thermal-centreline error and
is not the default:
<https://www.nist.gov/publications/fundamental-equations-state-parahydrogen-normal-hydrogen-and-orthohydrogen>,
<https://coolprop.org/fluid_properties/fluids/ParaHydrogen.html>.

**McCarty, R.D., Hord, J. and Roder, H.M. (1981).** *Selected Properties of
Hydrogen (Engineering Design Data)*, NBS Monograph 168. Its ortho/para section
reports that uncatalysed conversion is very slow, including a half-life over
one year at liquid-air temperatures. This rules out instantaneous spin
equilibrium across a millisecond jet and supports the fixed-composition
sensitivity in
[`prereg-hydrogen-spin-isomer-enthalpy.md`](prereg-hydrogen-spin-isomer-enthalpy.md):
<https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/nbsmonograph168.pdf>.

**Rattigan, W., Vizma, J. and Welch, J. (2026).** *D4.6 release into cold room
TCS tests*, ELVHYS grant 101101381, final submission 9 February 2026. Read with
the public ELVHYS WP4.2 dataset, DOI 10.18710/JXJP0H. The report documents the
1 m3 enclosure, test matrix and sensor layout, but neither source provides a
measured H2 mass-flow channel for dispersion Tests 10/11. The report and
archive sensor sheet conflict on horizontal nozzle elevation (200 versus
250 mm), so the screen stops before model scoring:
<https://cordis.europa.eu/project/id/101101381/results>,
<https://doi.org/10.18710/JXJP0H>.

**Jordan, T., Friedrich, A., Veser, A. and Kotchourko, N. (2019).**
*PRESLHY Experiment series E3.1 (Cryogenic Hydrogen Blow-down/Discharge)
results -- part A "high pressure"*. KITopen doi:10.5445/IR/1000096833;
RADAR4KIT archive doi:10.35097/1187, CC BY 4.0. The public archive contains
80 K/300 K high-pressure result files for 0.5/1/2/4 mm nozzles and a D=4 mm,
80 K camera/BOS image set. Its own description says the visible particles may
be nozzle ice or condensed hydrogen, so it supplies no N2/O2 particle-size,
phase-inventory, or gas-particle-slip closure. The local raw archive is ignored
and not redistributed; its intake/scope audit is
`public-preslhy-e31-intake.md`.

**Cirrone, D., Makarov, D., Kashkarov, S., Friedrich, A. and Molkov, V.
(2023).** "Physical model of non-adiabatic blowdown of cryo-compressed
hydrogen storage tanks." *International Journal of Hydrogen Energy*, 48(90),
35387--35406. doi:10.1016/j.ijhydene.2023.05.182. The open published version
uses the real-gas NIST/CoolProp properties, tank and discharge-pipe thermal
models, and an under-expanded nozzle construction; it reports validation
against sixteen PRESLHY tests (80--310 K, 0.6--20 MPa, 0.5--4 mm) with
test-specific Cd in 0.6--0.8. DEGALI's opt-in implementation is deliberately
smaller and documents its unresolved wall/withdrawal closures in
`cryocompressed-blowdown-status.md`.

**Cirrone, D., Makarov, D. V., Kuznetsov, M., Friedrich, A. and Molkov, V.
(2022).** "Effect of heat transfer through the release pipe on simulations of
cryogenic hydrogen jet fires and hazard distances." *International Journal of
Hydrogen Energy*, 47(50), 21596--21611. doi:10.1016/j.ijhydene.2022.04.276.
The open conference version reports geometry-dependent pipe warming and is
used only to justify a declared pipe-thermal sensitivity, never a universal
source-temperature correction.

**Shangguan, S., Wang, L., Shi, R., Li, Z., Xu, Z., Tan, H., Li, Y. and Lei,
G. (2025).** "Experimental investigation on
solidification characteristics of air-like nitrogen-oxygen mixtures in liquid
hydrogen." *Cryogenics*, article 104237.
doi:10.1016/j.cryogenics.2025.104237. The reported approximately 0.5--1.2 mm
particles were formed by injecting gas into bulk LH2. They do not provide a
nucleation-size closure for a freely expanding cryogenic jet and therefore do
not justify replacing the existing 1/10/100 um airborne-particle sensitivity.

**Steele, R.L. (1968).** *The thermodynamic properties of nitrogen-oxygen
mixtures*, MIT M.S./Naval Engineer thesis. A public 192-page archive copy was
read on 2026-09-17. It constructs N2/O2 liquid--vapour equilibrium
enthalpy/composition diagrams at 1, 5, 10 and 20 atm, but its atmospheric
coverage starts at 77.36 K. It is therefore provenance for rejecting a
pure-component enthalpy average and **not** a 59--63 K solid/liquid-mixture
closure. Local copy and any extracted numerical material remain ignored under
`reference/`; no external data are distributed:
<https://archive.org/download/thermodynamicpro00stee/thermodynamicpro00stee.pdf>.

**Hanna, S.R., Chang, J.C. and Strimaitis, D.G.** The acceptance criteria used
throughout — MG between 0.7 and 1.3, VG below 1.6, FAC2 above 0.5 — are theirs
as adopted by SMEDIS and the LNG model evaluation protocol. The primary paper
has not been read; the criteria are taken as they appear in those protocols.

**Winters and Houf**, whose plume model is in HyRAM; **Kamotani and Greber**,
whose trajectory correlation is in `SETJET`; **Pratte and Baines**, whose
development length it uses; **Birch and co-workers**, and **Yüceil and
Ötügen**, for the notional nozzle; **Ooms**, for the jet formulation
`JETPLU` implements; **Briggs**, for the lift-off thresholds. All are cited
through the reports that use them.

## The reason for this page

A parallel reimplementation of SLAB audited its own bibliography against the
primary sources and found a reference that does not exist — a report attributed
to a year in which it was not published — alongside four notation
inconsistencies and three citations that did not support the claim made from
them. The list had been complete in the sense that every entry was cited and
every citation listed; what was wrong was underneath that.
