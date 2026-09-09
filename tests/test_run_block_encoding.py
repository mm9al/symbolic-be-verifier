from evaluation import run_block_encoding


def test_block_encoding_runner_defaults_to_paper_threshold_without_residual_floor(monkeypatch):
    monkeypatch.setattr("sys.argv", ["run_block_encoding.py"])

    args = run_block_encoding._parse_args()

    assert args.check_mode == "paper"
    assert args.block_encoding_residual_tolerance == 0.0
