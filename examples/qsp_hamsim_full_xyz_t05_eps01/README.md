# Full Hamiltonian Simulation: H = (X + Y + Z)/3

This example verifies the degree-3 full Hamiltonian-simulation QSP block on the Hermitian base

```text
H = (X + Y + Z)/3
```

The QASM uses abstract `UH`/`UHdg` block-encoding gates. The concrete Hamiltonian is supplied to the symbolic verifier through `--base`.

Run from the repository root:

```bash
.venv/bin/python -m symbolic.verify examples/qsp_hamsim_full_xyz_t05_eps01/qsp_hamsim_full_xyz_t05_eps01_deg3.qasm \
  --ancillas 'q[0]' 'q[1]' 'q[2]' 'q[3]' \
  --systems 'q[4]' \
  --base '(X + Y + Z)/3' \
  --expected-polynomial '0.24991946353954456 - 0.12497982382931781*i*x - 0.030604023458682638*x^2 + 0.005127459989174488*i*x^3' \
  --hermitian-base \
  --result-only
```

Expected output:

```text
PASS
```
