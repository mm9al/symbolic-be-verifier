"""Sound truncation accounting, conservative acceptance, and mode isolation."""
import csv
from itertools import product
import json

import pytest
import sympy as sp

from evaluation.pennylane import bounded_pruning as sparse
from evaluation.pennylane import run_block_encoding as runner
from evaluation.pennylane.common import BENCHMARKS, sha256
from symbolic.branch_state import BranchState, Gate
from symbolic.expr import OpExpr
from symbolic.interval import acb, arb, ctx, to_acb
from symbolic.verify import verify_qasm_file


def artifact(n):
    directory = BENCHMARKS/'ising'/f'ising_periodic_n{n}'
    path = directory/'metadata.json'
    metadata = json.loads(path.read_text())
    keys = ('benchmark_id', 'model', 'graph', 'n_system', 'n_terms', 'target_lcu_terms',
            'n_ancilla', 'nominal_alpha', 'gate_count')
    row = dict(**{key: str(metadata[key]) for key in keys},
        qasm_path=str((directory/'block_encoding.qasm').relative_to(BENCHMARKS)),
        metadata_path=str(path.relative_to(BENCHMARKS)), metadata_sha256=sha256(path))
    return metadata, row, directory/'block_encoding.qasm'


def normalized_frobenius_difference(left, right):
    squared = arb(0)
    for bits in product((0, 1), repeat=left.num_ancillas):
        a, b = left.branch(bits).terms, right.branch(bits).terms
        for word in set(a) | set(b):
            squared += (to_acb(a.get(word, 0))-to_acb(b.get(word, 0))).abs_upper()**2
    return squared.sqrt().abs_upper()


def comparison_roundoff(left, right):
    # Two separately propagated enclosures are not correlated. Their difference
    # includes both radii even when the exact operators agree. Two times this
    # explicit radius sum bounds that comparison's excess over the exact error.
    radius = arb(0)
    for state in (left, right):
        for operator in state.branches.values():
            for coefficient in operator.terms.values():
                ball = to_acb(coefficient)
                radius += ball.real.rad()+ball.imag.rad()
    return (2*radius).upper()


def test_sum_over_all_branches_and_complex_coefficients_is_charged():
    with ctx.workdps(100):
        tiny = acb(arb('3e-30'), arb('4e-30'))
        state = BranchState(1, 1, {(0,): OpExpr({('I',): 1, ('X',): tiny}),
                                  (1,): OpExpr({('Y',): -tiny})})
        ledger = sparse.ErrorLedger('1e-24', '1e-14')
        result = ledger.prune(state, 1)
        assert result.branch((0,)).terms == {('I',): 1}
        assert result.branch((1,)).is_zero()
        assert ledger.dropped_terms == 2
        assert ledger.error >= arb('1e-29').upper()
        assert ledger.error < arb('1.01e-29')
        assert normalized_frobenius_difference(state, result) <= ledger.error


def test_discarded_nonzero_ancilla_branch_can_return_to_selected_branch():
    # A ledger looking only at the currently selected branch would charge zero.
    # X on the ancilla later moves the entire discarded operator into that branch.
    with ctx.workdps(100):
        tiny = sp.Rational(1, 2**40)
        # This column is an isometry and can be reached by a joint unitary.
        leading = sp.sqrt(1-tiny**2)
        state = BranchState(1, 1, {(0,): OpExpr({('I',): leading}), (1,): OpExpr({('Z',): tiny})})
        ledger = sparse.ErrorLedger('1e-10', '1e-8')
        pruned = ledger.prune(state, 1)
        assert (state.top_left()-pruned.top_left()).is_zero()
        assert ledger.error >= to_acb(tiny).abs_upper()
        original_after = state.apply(Gate('x', (0,)), ancillas=[0], systems=[1])
        pruned_after = pruned.apply(Gate('x', (0,)), ancillas=[0], systems=[1])
        difference = original_after.top_left()-pruned_after.top_left()
        assert not difference.is_zero()
        # Z is unitary, so the operator norm of tiny*Z equals |tiny| exactly.
        assert to_acb(difference.terms[('Z',)]).abs_upper() <= ledger.error


def test_zero_containing_balls_are_not_automatically_discarded():
    state = BranchState(1, 1, {(0,): OpExpr({('I',): 1, ('X',): acb(arb('0 +/- 1e-3'))})})
    ledger = sparse.ErrorLedger('1e-24', '1e-14')
    assert ledger.prune(state, 1) is state
    assert ledger.error.is_zero()


def test_budget_exhaustion_keeps_entire_batch_and_zero_budget_disables_pruning():
    state = BranchState(1, 1, {(0,): OpExpr({('I',): 1, ('X',): sp.Rational(1, 10)})})
    for budget in ('0', '0.05'):
        ledger = sparse.ErrorLedger('0.2', budget)
        assert ledger.prune(state, 1) is state
        assert ledger.error.is_zero()
        assert ledger.dropped_terms == 0
    assert ledger.skipped_batches == 1


def test_later_batches_use_the_remaining_budget():
    state = BranchState(1, 1, {(0,): OpExpr({('I',): 1, ('X',): sp.Rational(1, 10)})})
    ledger = sparse.ErrorLedger('0.2', '0.15')
    assert ('X',) not in ledger.prune(state, 1).top_left().terms
    assert ledger.prune(state, 2) is state
    assert ledger.dropped_terms == 1
    assert ledger.skipped_batches == 1
    assert ledger.error <= ledger.budget_ball.lower()


def test_recentering_charges_full_radii_and_keeps_exact_dyadic_centers():
    with ctx.workdps(100):
        coefficient = acb(arb('0.5 +/- 1e-8'), arb('0.25 +/- 2e-8'))
        state = BranchState(1, 1, {(0,): OpExpr({('I',): coefficient}),
                                  (1,): OpExpr({('Z',): coefficient})})
        ledger = sparse.ErrorLedger('1e-24', '1e-6')
        centered = ledger.recenter(state, 1)
        assert ledger.dropped_terms == 0
        assert ledger.recentered_coefficients == 2
        assert ledger.error >= arb('6e-8').upper()
        assert ledger.roundoff_error == ledger.error
        for operator in centered.branches.values():
            for value in operator.terms.values():
                assert value.real.is_exact() and value.imag.is_exact()
        assert normalized_frobenius_difference(state, centered) <= ledger.error
        capped = sparse.ErrorLedger('1e-24', '1e-9')
        assert capped.recenter(state, 1) is state
        assert capped.error.is_zero()
        assert capped.skipped_recentering_batches == 1


@pytest.mark.parametrize('real_radius,imag_radius', [(0, 2**-20), (2**-20, 0),
                                                   (2**-20, 2**-21)])
def test_recentering_bounds_complex_corners_in_both_coordinates(real_radius, imag_radius):
    with ctx.workdps(100):
        coefficient = acb(arb('0.375', real_radius), arb('-0.25', imag_radius))
        state = BranchState(1, 1, {(1,): OpExpr({('X',): coefficient})})
        ledger = sparse.ErrorLedger('1e-24', '1e-3')
        centered = ledger.recenter(state, 1)
        midpoint = centered.branch((1,)).terms[('X',)]
        radius_sum = (coefficient.real.rad()+coefficient.imag.rad()).upper()
        assert ledger.error >= radius_sum
        for real_sign, imag_sign in product((-1, 1), repeat=2):
            corner = acb(coefficient.real.mid()+real_sign*coefficient.real.rad(),
                         coefficient.imag.mid()+imag_sign*coefficient.imag.rad())
            assert (corner-midpoint).abs_upper() <= ledger.error


def test_recentering_stabilizes_a_long_exact_inverse_circuit(tmp_path):
    # Each forward/inverse pair is exactly I for the serialized decimal angle.
    path = tmp_path/'inverse_pairs.qasm'
    forward = ['rz(1.23456789) q[1];', 'rx(0.91827364) q[1];']
    reverse = ['rx(-0.91827364) q[1];', 'rz(-1.23456789) q[1];']
    path.write_text('OPENQASM 2.0;\nqreg q[2];\n'+'\n'.join((forward+reverse)*60))
    with ctx.workdps(30):
        result, _, report = sparse.verify_with_pruning(path, ancillas=[0], systems=[1], precision=30)
        difference = result.final_state.top_left()-OpExpr.identity()
        norm = sum((to_acb(c).abs_upper() for c in difference.terms.values()), arb(0))
        eta = arb(report['operator_error_upper']).abs_upper()
        assert norm <= eta
        assert eta < arb('1e-20')
        assert report['recentered_coefficients'] > 100


@pytest.mark.parametrize('threshold,budget', [('0', '1'), ('-1', '1'), ('1', '-1')])
def test_invalid_configuration(threshold, budget):
    with pytest.raises(ValueError):
        sparse.ErrorLedger(threshold, budget)


@pytest.mark.parametrize('statement', ['cx q[0],q[0];', 'h q[0],q[1];',
                                     'rx q[1];', 'h(1) q[1];', 'uh q[0],q[1];'])
def test_invalid_or_opaque_operations_cannot_use_the_unitary_error_proof(tmp_path, statement):
    path = tmp_path/'invalid.qasm'
    path.write_text('OPENQASM 2.0;\nqreg q[2];\n'+statement)
    with pytest.raises(ValueError):
        sparse.verify_with_pruning(path, ancillas=[0], systems=[1], precision=100)


def test_overlapping_registers_cannot_use_the_unitary_error_proof(tmp_path):
    path = tmp_path/'overlap.qasm'
    path.write_text('OPENQASM 2.0;\nqreg q[2];\nx q[0];')
    with pytest.raises(ValueError, match='disjoint registers'):
        sparse.verify_with_pruning(path, ancillas=[0], systems=[0], precision=100)


@pytest.mark.parametrize('n', [2, 3, 4])
def test_native_all_branch_errors_obey_the_global_ledger(n):
    metadata, _, path = artifact(n)
    original_hash = sha256(path)
    with ctx.workdps(100):
        result, _, report = sparse.verify_with_pruning(path, ancillas=metadata['ancillas'],
            systems=metadata['systems'], precision=100)
        reference = verify_qasm_file(path, ancillas=metadata['ancillas'], systems=metadata['systems'],
                                     qsp_calculation_precision=100)
        eta = arb(report['operator_error_upper'])
        difference = normalized_frobenius_difference(reference.final_state, result.final_state)
        assert difference <= (eta.abs_upper()+comparison_roundoff(reference.final_state, result.final_state)).upper()
        assert eta <= arb(sparse.DEFAULT_BUDGET)
        assert sha256(path) == original_hash


def test_error_bound_survives_noncommuting_unitary_suffix_and_ancilla_mixing(tmp_path):
    path = tmp_path/'noncommuting.qasm'
    path.write_text('''OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
rx(1e-28) q[1];
ry(0.71) q[0];
cy q[0],q[2];
rz(1.37) q[1];
h q[0];
ry(-0.43) q[2];
cx q[1],q[0];
''')
    with ctx.workdps(100):
        result, _, report = sparse.verify_with_pruning(path, ancillas=[0], systems=[1, 2], precision=100)
        reference = verify_qasm_file(path, ancillas=[0], systems=[1, 2], qsp_calculation_precision=100)
        eta = arb(report['operator_error_upper'])
        assert eta > 0
        allowance = comparison_roundoff(reference.final_state, result.final_state)
        assert normalized_frobenius_difference(reference.final_state, result.final_state) <= (eta+allowance).upper()


def test_scaled_error_is_charged_to_original_operator_acceptance():
    # The sparse residual alone passes, but the original circuit is uncertified.
    certificate = sparse.residual_certificate(arb('1e-10'), arb(10), arb(4), '1e-8', arb('2e-9'))
    assert not certificate['passed']
    assert certificate['scaled_truncation_error_upper'] >= arb('2e-8')
    assert sparse.residual_certificate(arb('1e-10'), arb(10), arb(4), '1e-8', arb('1e-12'))['passed']
    assert not sparse.residual_certificate(arb(0), arb(-1), arb(4), '1e-8', arb(0))['passed']


def test_pruning_cannot_turn_unknown_exact_hermiticity_into_a_proof():
    actual = OpExpr({('I',): 1})
    proof = {'status': 'UNPROVEN'}
    result = sparse.classify_with_pruning(actual, proof, '0', arb('1e-13'))
    assert result['status'] == result['exact_status'] == 'UNPROVEN'
    assert sparse.classify_with_pruning(actual, proof, '1e-12', arb('1e-13'))['status'] == 'PASS_WITHIN_DELTA'
    # Covers missing/discarded Pauli coefficients, including an entirely empty block.
    assert sparse.classify_with_pruning(OpExpr.zero(), proof, '1e-12', arb('2e-12'))['status'] == 'FAIL_DELTA'


def test_hermiticity_strict_boundary_and_nonhermitian_witness_include_error():
    proof = {'status': 'UNPROVEN'}
    actual = OpExpr({('X',): acb(0, arb('0.75'))})
    assert sparse.classify_with_pruning(actual, proof, '1', arb('0.25'))['status'] == 'FAIL_DELTA'
    result = sparse.classify_with_pruning(actual, proof, '1', arb('0.125'))
    assert result['status'] == 'PASS_WITHIN_DELTA'
    assert result['exact_status'] == 'FAIL_NON_HERMITIAN'
    assert sparse.classify_with_pruning(actual, proof, '0', arb(1))['exact_status'] == 'UNPROVEN'


@pytest.mark.parametrize('delta,status', [('1e-12', 'ACCEPT'), ('0', 'REJECT')])
def test_end_to_end_keeps_decimal_circuit_and_correct_acceptance(tmp_path, delta, status):
    _, row, path = artifact(3)
    result = runner.evaluate(BENCHMARKS, row, tmp_path, '1e-8', 100, delta, bounded_pruning=True)
    assert result['status'] == status, result
    assert result['exact_hermiticity_status'] == 'FAIL_NON_HERMITIAN'
    report = json.loads((tmp_path/result['report_path']).read_text())
    assert report['qasm_sha256'] == sha256(path)
    assert report['bounded_pruning']['original_qasm_verified']
    assert 'alpha*eta' in report['acceptance_rule']


def test_cli_has_distinct_default_directory_and_extra_schema(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, 'HERE', tmp_path)
    monkeypatch.setattr(runner, 'ROOT', tmp_path)
    for name in ('run_block_encoding.py', 'common.py', 'local_summaries.py', 'bounded_pruning.py'):
        (tmp_path/name).write_text('# dispatch fixture\n')
    monkeypatch.setattr('sys.argv', ['run_block_encoding.py', '--manifest',
                                   str(BENCHMARKS/'delta_calibration_manifest.csv'),
                                   '--sizes', '2', '--bounded-pruning',
                                   '--timeout-sec', '0'])
    # Keep this a dispatch/schema test; real execution is tested above.
    def fake(*args):
        assert args[8] is True
        return dict(benchmark_id='ising_periodic_n2', status='ACCEPT', success=True)
    monkeypatch.setattr(runner, 'run_case', fake)
    assert runner.main() == 0
    output = tmp_path/'results_bounded_pruning'
    with (output/'block_encoding_results.csv').open() as handle:
        assert csv.DictReader(handle).fieldnames == runner.FIELDS+runner.SUMMARY_FIELDS+runner.PRUNING_FIELDS
    assert not (tmp_path/'results').exists()
    assert not (tmp_path/'results_local_summaries').exists()
