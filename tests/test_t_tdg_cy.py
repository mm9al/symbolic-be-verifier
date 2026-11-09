"""Check branch updates against independent exact computational-basis matrices."""

from itertools import permutations, product

import pytest
import sympy as sp

from symbolic.branch_state import BranchState, Gate, UnsupportedGateError
from symbolic.expr import OpExpr
from symbolic.verify import PASS_EXACT, verify_qasm_file
from symbolic.word_expr import WordExpr


PAULIS = {
    "I": sp.eye(2),
    "X": sp.Matrix([[0, 1], [1, 0]]),
    "Y": sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    "Z": sp.diag(1, -1),
}
LOCAL = {
    "t": sp.diag(1, sp.exp(sp.I * sp.pi / 4)),
    "tdg": sp.diag(1, sp.exp(-sp.I * sp.pi / 4)),
    "cy": sp.diag(sp.eye(2), PAULIS["Y"]),
}
ANCILLAS = (7, 2)
SYSTEMS = (9, 4)
WIRE_ORDER = ANCILLAS + SYSTEMS


def dense_operator(operator):
    result = sp.zeros(2**operator.num_qubits)
    for word, coefficient in operator.terms.items():
        result += coefficient * sp.kronecker_product(*(PAULIS[axis] for axis in word))
    return result


def dense_branches(state):
    return sp.Matrix.vstack(
        *(dense_operator(state.branch(bits)) for bits in product((0, 1), repeat=state.num_ancillas))
    )


def embedded_matrix(local, operands, wire_order=WIRE_ORDER):
    """Place a literal gate matrix on specified wires, using basis columns."""
    matrix = sp.zeros(2 ** len(wire_order))
    positions = [wire_order.index(wire) for wire in operands]
    for column, bits in enumerate(product((0, 1), repeat=len(wire_order))):
        local_column = int("".join(str(bits[position]) for position in positions), 2)
        for local_row, changed_bits in enumerate(product((0, 1), repeat=len(operands))):
            output = list(bits)
            for position, bit in zip(positions, changed_bits):
                output[position] = bit
            row = int("".join(map(str, output)), 2)
            matrix[row, column] = local[local_row, local_column]
    return matrix


def sample_state(sparse):
    # Independent, noncommuting operators on every branch expose right-multiply
    # errors, control/target swaps, and the relative phases between branches.
    branches = {
        bits: OpExpr(
            {
                ("I", "I"): index + 1,
                ("X", "Z"): index + 2 * sp.I,
                ("Y", "X"): 2 - index,
                ("Z", "Y"): sp.Rational(index + 1, 3),
            },
            num_qubits=2,
        )
        for index, bits in enumerate(product((0, 1), repeat=2))
        if not sparse or bits == (1, 0)
    }
    return BranchState(2, 2, branches)


GATES = (
    [Gate(name, (wire,)) for name in ("t", "tdg") for wire in WIRE_ORDER]
    + [Gate("cy", operands) for operands in permutations(WIRE_ORDER, 2)]
)


@pytest.mark.parametrize("sparse", (False, True))
@pytest.mark.parametrize("gate", GATES, ids=str)
def test_branch_update_matches_exact_dense_matrix(gate, sparse):
    before = sample_state(sparse)
    after = before.apply(gate, ancillas=ANCILLAS, systems=SYSTEMS)
    expected = embedded_matrix(LOCAL[gate.name], gate.qubits) * dense_branches(before)
    difference = dense_branches(after) - expected

    assert all(sp.simplify(sp.expand_complex(entry)) == 0 for entry in difference)


@pytest.mark.parametrize("wire", WIRE_ORDER)
def test_t_and_tdg_are_exact_inverses_and_t_squared_is_s(wire):
    before = sample_state(False)
    t = before.apply(Gate("t", (wire,)), ancillas=ANCILLAS, systems=SYSTEMS)
    inverse = t.apply(Gate("tdg", (wire,)), ancillas=ANCILLAS, systems=SYSTEMS)
    squared = t.apply(Gate("t", (wire,)), ancillas=ANCILLAS, systems=SYSTEMS)
    s = before.apply(Gate("s", (wire,)), ancillas=ANCILLAS, systems=SYSTEMS)

    for bits in product((0, 1), repeat=2):
        assert inverse.branch(bits).equals(before.branch(bits))
        assert squared.branch(bits).equals(s.branch(bits))


@pytest.mark.parametrize("name,phase", [("t", (1 + sp.I) / sp.sqrt(2)), ("tdg", (1 - sp.I) / sp.sqrt(2))])
def test_t_phase_update_in_word_mode(name, phase):
    before = BranchState(1, 1, {(0,): WordExpr({("H",): 2}), (1,): WordExpr({("G",): 3})})
    after = before.apply(Gate(name, (0,)))

    assert after.b0.equals(before.b0)
    assert after.b1.equals(before.b1.scale(phase))


def test_cy_ancilla_branch_phases_in_word_mode():
    before = BranchState(2, 1, {(1, 0): WordExpr({("H",): 2}), (1, 1): WordExpr({("G",): 3})})
    after = before.apply(Gate("cy", (0, 1)), ancillas=(0, 1), systems=(2,))

    assert after.branch((1, 1)).equals(before.branch((1, 0)).scale(sp.I))
    assert after.branch((1, 0)).equals(before.branch((1, 1)).scale(-sp.I))


@pytest.mark.parametrize("gate", [Gate("t", (4,)), Gate("tdg", (4,)), Gate("cy", (0, 4)), Gate("cy", (4, 0)), Gate("cy", (4, 5))])
def test_explicit_system_gates_remain_rejected_in_word_mode(gate):
    with pytest.raises(UnsupportedGateError, match="word mode"):
        BranchState.initial(2, 1, expression_kind="word").apply(gate, ancillas=(0,), systems=(4, 5))


@pytest.mark.parametrize("operands", [(0, 0), (1, 1), (0, 8), (8, 0), (8, 1), (1, 8)])
def test_cy_rejects_duplicate_or_unknown_operands(operands):
    with pytest.raises(UnsupportedGateError, match="distinct|known"):
        BranchState.initial().apply(Gate("cy", operands))


@pytest.mark.parametrize("gate", [Gate("t", (0,), parameter=sp.pi), Gate("tdg", (1,), parameter=sp.pi), Gate("cy", (0, 1), parameter=sp.pi)])
def test_fixed_gates_reject_parameters(gate):
    with pytest.raises(UnsupportedGateError, match="Parameterized gate"):
        BranchState.initial().apply(gate)



def test_qasm_t_and_tdg_preserve_absolute_phase(tmp_path):
    path = tmp_path / "phase.qasm"
    path.write_text('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[2];\nh q[0];\nt q[0];\nh q[0];\ntdg q[1];\n')
    result = verify_qasm_file(path, expected="((1+exp(i*pi/4))/2)*((1+exp(-i*pi/4))*I/2+(1-exp(-i*pi/4))*Z/2)")

    assert result.status == PASS_EXACT
