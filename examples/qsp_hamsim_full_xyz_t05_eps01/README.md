# Full Hamiltonian Simulation: H = (X + Y + Z)/3

This example verifies the degree-3 full Hamiltonian-simulation QSP block on the Hermitian base

```text
H = (X + Y + Z)/3
```

The QASM uses abstract `UH`/`UHdg` block-encoding gates. The verifier first
extracts the symbolic polynomial, then checks it against the Hamiltonian
simulation target on `[-1, 1]`.

Run from the repository root:

```bash
.venv/bin/python -m symbolic.verify examples/qsp_hamsim_full_xyz_t05_eps01/qsp_hamsim_full_xyz_t05_eps01_deg3.qasm \
  --ancillas 'q[0]' 'q[1]' 'q[2]' 'q[3]' \
  --systems 'q[4]' \
  --hamsim-tau 0.5 \
  --hamsim-epsilon 0.1 \
  --hamsim-scale 0.25 \
  --hermitian-base \
  --result-only
```

Expected output:

```text
PASS
```
