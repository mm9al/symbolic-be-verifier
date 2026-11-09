"""Exact equivalence, fallback behavior, branch updates, and opt-in isolation."""
from collections import Counter
import csv
from itertools import permutations, product
import json

import pytest
import sympy as sp

from evaluation.pennylane import local_summaries as local
from evaluation.pennylane import run_block_encoding as runner
from evaluation.pennylane.common import BENCHMARKS, sha256
from symbolic.branch_state import BranchState, Gate
from symbolic.expr import OpExpr
from symbolic.interval import arb, ctx, to_acb
from symbolic.qasm_parser import parse_qasm_file
from symbolic.verify import verify_qasm_file


def artifact(n):
    directory = BENCHMARKS/"ising"/f"ising_periodic_n{n}"
    path = directory/"metadata.json"
    metadata = json.loads(path.read_text())
    keys = ("benchmark_id", "model", "graph", "n_system", "n_terms", "target_lcu_terms",
            "n_ancilla", "nominal_alpha", "gate_count")
    row = dict(**{key: str(metadata[key]) for key in keys},
        qasm_path=str((directory/"block_encoding.qasm").relative_to(BENCHMARKS)),
        metadata_path=str(path.relative_to(BENCHMARKS)), metadata_sha256=sha256(path))
    return metadata, row, directory/"block_encoding.qasm"


def instantiate(template, wires=(7, 2, 11)):
    return [Gate(name, tuple(wires[role] for role in roles)) for name, roles in template.operations]


@pytest.mark.parametrize("template", local.TEMPLATES, ids=lambda t: t.name)
def test_templates_have_exact_all_column_certificates(template):
    proof = local.certify_template(template)
    assert proof["status"] == "PASS_EXACT"
    assert proof["input_columns"] == 8
    assert proof["absolute_phase_checked"]


@pytest.mark.parametrize("template", local.TEMPLATES, ids=lambda t: t.name)
def test_changed_t_phase_fails_the_exact_certificate(template):
    operations = list(template.operations)
    index = next(i for i, (name, _) in enumerate(operations) if name == "tdg")
    operations[index] = ("t", operations[index][1])
    corrupt = local.Template("corrupt", template.primitive, tuple(operations))
    with pytest.raises(ValueError, match="Exact equivalence proof failed"):
        local.certify_template(corrupt)


def test_global_minus_sign_is_not_accepted_as_equivalent():
    # (X Z)^2 = -I; synthesize X as H T^4 H and Z as T^4.
    z = (("t", (0,)),)*4
    x = (("h", (0,)),)+z+(("h", (0,)),)
    negative_identity = z+x+z+x
    template = local.Template("wrong_absolute_phase", "ccx", local.PENNYLANE_CCX+negative_identity)
    with pytest.raises(ValueError, match="Exact equivalence proof failed"):
        local.certify_template(template)


@pytest.mark.parametrize('power', [1, 2])
def test_global_pi_over_four_and_i_phases_are_not_accepted_as_equivalent(power):
    # T^k X T^k X = exp(i*k*pi/4) I, with X = H T^4 H.
    t = (("t", (0,)),)
    x = (("h", (0,)),)+t*4+(("h", (0,)),)
    scalar_phase = t*power+x+t*power+x
    template = local.Template("wrong_complex_global_phase", "ccx", local.PENNYLANE_CCX+scalar_phase)
    with pytest.raises(ValueError, match="Exact equivalence proof failed"):
        local.certify_template(template)


@pytest.mark.parametrize("template", local.TEMPLATES, ids=lambda t: t.name)
def test_match_supports_remapped_wires_and_interleaved_exact_identities(template):
    original = instantiate(template)
    gates = [item for gate in original for item in (gate, Gate("rz", (11,), sp.Integer(0)))]
    plan = local.compile_summaries(gates, ancillas=(2, 7), systems=(11,))
    assert plan.report["summary_count"] == 1
    assert plan.report["removed_identity_gates"] == len(original)
    assert len(plan.steps) == 1
    assert plan.steps[0].qubits == (7, 2, 11)


def test_near_zero_rotation_breaks_the_pattern_and_remains_in_the_plan():
    gates = instantiate(local.TEMPLATES[0])
    tiny = Gate("rz", (11,), arb("1e-30"), raw="rz(1e-30) q[11];")
    gates.insert(3, tiny)
    plan = local.compile_summaries(gates, ancillas=(7, 2), systems=(11,))
    assert plan.report["summary_count"] == 0
    assert plan.report["removed_identity_gates"] == 0
    assert tiny in plan.steps


@pytest.mark.parametrize("gate", [Gate("rz", (99,), sp.Integer(0)),
                                  Gate("rz", (7, 11), sp.Integer(0))])
def test_identity_removal_does_not_hide_invalid_gate_operands(gate):
    plan = local.compile_summaries([gate], ancillas=(7, 2), systems=(11,))
    assert plan.steps == (gate,)
    with pytest.raises(ValueError):
        local.execute_plan(plan, ancillas=(7, 2), systems=(11,))


def test_unsupported_system_control_falls_back_to_original_gates():
    gates = instantiate(local.TEMPLATES[0], (11, 2, 7))
    plan = local.compile_summaries(gates, ancillas=(7, 2), systems=(11,))
    assert plan.report["summary_count"] == 0
    assert list(plan.steps) == gates


PAULIS = {"I": sp.eye(2), "X": sp.Matrix([[0, 1], [1, 0]]),
          "Y": sp.Matrix([[0, -sp.I], [sp.I, 0]]), "Z": sp.diag(1, -1)}


def dense_operator(operator):
    result = sp.zeros(2**operator.num_qubits)
    for word, coefficient in operator.terms.items():
        result += coefficient*sp.kronecker_product(*(PAULIS[p] for p in word))
    return result


def dense_branches(state):
    return sp.Matrix.vstack(*(dense_operator(state.branch(bits))
                              for bits in product((0, 1), repeat=state.num_ancillas)))


def literal_controlled_matrix(name, operands, wire_order):
    # Independent literal basis permutation/phase, without branch update helpers.
    positions = [wire_order.index(wire) for wire in operands]
    matrix = sp.zeros(2**len(wire_order))
    for column, bits in enumerate(product((0, 1), repeat=len(wire_order))):
        output = list(bits)
        phase = 1
        if bits[positions[0]] and bits[positions[1]]:
            if name == "ccx":
                output[positions[2]] ^= 1
            elif bits[positions[2]]:
                phase = -1
        row = int("".join(map(str, output)), 2)
        matrix[row, column] = phase
    return matrix


@pytest.mark.parametrize("template", local.TEMPLATES, ids=lambda t: t.name)
@pytest.mark.parametrize("wires", list(permutations((7, 2, 11))))
def test_matching_and_summary_embedding_preserve_all_role_permutations(template, wires):
    ancillas, systems = (11, 7, 2), (19,)
    plan = local.compile_summaries(instantiate(template, wires), ancillas=ancillas, systems=systems)
    assert plan.report["summary_count"] == 1
    assert plan.steps[0].qubits == wires
    state = BranchState(3, 1, {bits: OpExpr({("I",): i+1, ("X",): i+sp.I})
                                for i, bits in enumerate(product((0, 1), repeat=3))})
    actual = local.apply_summary(state, plan.steps[0], ancillas=ancillas, systems=systems)
    expected = literal_controlled_matrix(template.primitive, wires, ancillas+systems)*dense_branches(state)
    assert all(sp.expand(entry) == 0 for entry in dense_branches(actual)-expected)


@pytest.mark.parametrize("name", ["ccx", "ccz"])
@pytest.mark.parametrize("target", [10, 11])
def test_summary_updates_match_independent_matrices_on_all_branches(name, target):
    ancillas, systems = (7, 2, 10), (11, 4)
    state = BranchState(3, 2, {bits: OpExpr({("I", "I"): i+1, ("X", "Z"): i+sp.I,
                                          ("Z", "Y"): 2-i}, num_qubits=2)
                                for i, bits in enumerate(product((0, 1), repeat=3))})
    gate = local.SummaryGate(name, (7, 2, target), "test", 1, 15)
    actual = local.apply_summary(state, gate, ancillas=ancillas, systems=systems)
    expected = literal_controlled_matrix(name, gate.qubits, ancillas+systems)*dense_branches(state)
    assert all(sp.expand(entry) == 0 for entry in dense_branches(actual)-expected)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_native_whole_circuit_branch_enclosures_agree(n):
    metadata, _, path = artifact(n)
    original_hash = sha256(path)
    with ctx.workdps(100):
        optimized, report = local.verify_with_summaries(path, ancillas=metadata["ancillas"],
            systems=metadata["systems"], precision=100)
        reference = verify_qasm_file(path, ancillas=metadata["ancillas"], systems=metadata["systems"],
                                     qsp_calculation_precision=100)
        for key in set(optimized.final_state.branches) | set(reference.final_state.branches):
            left = optimized.final_state.branch(key)
            right = reference.final_state.branch(key)
            for word in set(left.terms) | set(right.terms):
                assert (to_acb(left.terms.get(word, 0))-to_acb(right.terms.get(word, 0))).contains(0)
        assert report["summary_count"] > 0
        assert report["approximation_error"] == "0"
        assert report["source_gate_count"] == metadata["gate_count"]
        assert sha256(path) == original_hash


def test_every_nonzero_decimal_rotation_is_retained():
    metadata, _, path = artifact(4)
    gates = parse_qasm_file(path, numeric_precision=100)
    plan = local.compile_summaries(gates, ancillas=metadata["ancillas"], systems=metadata["systems"])
    original = Counter(g.raw for g in gates if g.parameter is not None and not local.is_exact_zero(g.parameter))
    retained = Counter(g.raw for g in plan.steps if isinstance(g, Gate)
                       and g.parameter is not None and not local.is_exact_zero(g.parameter))
    assert retained == original


@pytest.mark.parametrize("delta,status", [("1e-12", "ACCEPT"), ("0", "REJECT")])
def test_summary_mode_retains_actual_decimal_nonhermiticity(tmp_path, delta, status):
    _, row, _ = artifact(3)
    result = runner.evaluate(BENCHMARKS, row, tmp_path, "1e-8", 100, delta, local_summaries=True)
    assert result["status"] == status, result
    assert result["exact_hermiticity_status"] == "FAIL_NON_HERMITIAN"
    report = json.loads((tmp_path/result["report_path"]).read_text())
    assert report["local_summaries"]["summary_count"] == 13
    assert report["qasm_sha256"] == sha256(BENCHMARKS/row["qasm_path"])


@pytest.mark.parametrize("enabled", [False, True])
def test_cli_keeps_original_csv_schema_unless_enabled(tmp_path, monkeypatch, enabled):
    arguments = ["run_block_encoding.py", "--manifest", str(BENCHMARKS/"delta_calibration_manifest.csv"),
                 "--sizes", "2", "--output-dir", str(tmp_path), "--timeout-sec", "0"]
    if enabled:
        arguments.append("--local-summaries")
    monkeypatch.setattr("sys.argv", arguments)
    assert runner.main() == 0
    with (tmp_path/"block_encoding_results.csv").open() as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == runner.FIELDS+(runner.SUMMARY_FIELDS if enabled else [])
        rows = list(reader)
    assert rows[0]["status"] == "ACCEPT"
    report = json.loads((tmp_path/rows[0]["report_path"]).read_text())
    assert ("local_summaries" in report) == enabled
    metadata = json.loads((tmp_path/"run_metadata.json").read_text())
    assert metadata["source_snapshots"].keys() == metadata["sources_sha256"].keys()
    for name, relative in metadata["source_snapshots"].items():
        assert sha256(tmp_path/relative) == metadata["sources_sha256"][name]
    manifest = tmp_path/metadata["manifest_snapshot"]
    assert manifest.read_bytes() == (BENCHMARKS/"delta_calibration_manifest.csv").read_bytes()
    assert sha256(manifest) == metadata["manifest_sha256"]
