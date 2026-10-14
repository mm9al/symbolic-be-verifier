from evaluation import run_block_encoding


def test_block_encoding_runner_has_only_the_paper_epsilon(monkeypatch):
    monkeypatch.setattr("sys.argv", ["run_block_encoding.py"])

    args = run_block_encoding._parse_args()

    assert args.check_mode == "paper"
    assert args.block_encoding_epsilon == 1e-8
    assert not hasattr(args, "block_encoding_residual_tolerance")

    fieldnames = run_block_encoding._fieldnames()
    assert "residual_threshold" in fieldnames
    assert "residual_numerical_tolerance" not in fieldnames
    assert "residual_acceptance_threshold" not in fieldnames
