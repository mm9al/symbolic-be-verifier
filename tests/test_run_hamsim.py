from evaluation import run_hamsim


def test_hamsim_runner_defaults_to_paper_approximation_on_fixed_epsilon_axis(monkeypatch):
    monkeypatch.setattr("sys.argv", ["run_hamsim.py"])

    args = run_hamsim._parse_args()

    assert args.check_mode == "approximation"
    assert args.axes == ["vary_t"]
