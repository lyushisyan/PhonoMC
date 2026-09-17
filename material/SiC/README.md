# 3C-SiC transport material

This material uses the **SiC-specific NEP interatomic potential**, with 5×5×5 force-constant supercells and a 31×31×31 q mesh. It is not a direct DFT force-constant dataset. The files match the material used for the Si/SiC interface and FinFET calculations.

- `POSCAR`: primitive cell consistent with the phonon table.
- `kappa-fbz.hdf5`: 29,791 full-zone q points, 6 branches, unit weights; temperatures 100–1000 K in 50 K increments.
- Scattering: `gamma`, `gamma_N`, `gamma_U`, and independent `gamma_isotope`. N+U reproduces `gamma`; isotope scattering is added separately to RTA and the Callaway resistive channel.
- Linewidth units: THz HWHM, with inverse lifetime `4*pi*gamma` in ps⁻¹. Frequencies are in THz and velocities use the phono3py THz·Å convention.
- The original 1 m boundary mean-free-path contribution is omitted from transport rates because the solver treats device boundaries explicitly.
- `validation.json`: archived full-zone export checks and bulk-conductivity comparisons. `provenance.json`: file hashes and dataset source identifiers.

The transport table uses a 0.01 THz frequency cutoff; unmodified phono3py frequencies are retained as `frequency_raw_phono3py`. See the validation record for the corresponding conductivity comparison. The source metadata inside the HDF5 is retained unchanged.

RTA supports local linewidth lookup within the supplied temperature range. Callaway additionally applies its own mode-grid and reference-temperature checks; having N/U arrays alone does not establish convergence for a particular geometry. Choose temperatures and boundary discretizations accordingly.
