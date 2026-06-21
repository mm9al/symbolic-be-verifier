#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ -z "${PYTHON:-}" && -x .venv/bin/python ]]; then
  PYTHON=.venv/bin/python
else
  PYTHON=${PYTHON:-python3}
fi

"$PYTHON" -m symbolic.verify examples/lcu_x_plus_z.qasm --expected "(X + Z)/2"

"$PYTHON" -m symbolic.verify examples/two_ancilla_hh.qasm \
  --ancillas 'q[0]' 'q[1]' \
  --systems 'q[2]' \
  --expected "I/2"

"$PYTHON" -m symbolic.verify examples/lcu_xx_plus_zz.qasm \
  --systems 1 2 \
  --expected "(X0 X1 + Z0 Z1)/2"

"$PYTHON" -m symbolic.verify examples/rz_ancilla_theta_sandwich.qasm \
  --expected "cos(theta/2) I"

echo "All smoke tests passed."
