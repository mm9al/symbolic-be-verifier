import sympy as sp
from flint import acb, arb, ctx

from symbolic.interval import to_acb, to_arb
from symbolic.scalar import cos_half, exp_plus_i_half, parse_scalar, scalar_simplify, sin_half


def test_parse_scalar_supported_angle_forms():
    theta = parse_scalar("theta")

    assert parse_scalar("pi") == sp.pi
    assert parse_scalar("pi/2") == sp.pi / 2
    assert parse_scalar("0.5*pi") == sp.pi / 2
    assert parse_scalar("-theta") == -theta
    assert parse_scalar("2*pi") == 2 * sp.pi


def test_scalar_simplification_examples():
    theta = sp.Symbol("theta")

    assert scalar_simplify(sp.cos(sp.pi / 4)) == sp.sqrt(2) / 2
    assert scalar_simplify(sp.sin(sp.pi / 4)) == sp.sqrt(2) / 2
    assert scalar_simplify(sp.exp(sp.I * sp.pi)) == -1
    assert scalar_simplify(sp.sin(theta / 2) ** 2 + sp.cos(theta / 2) ** 2) == 1


def test_rotation_exp_helpers_keep_symbolic_theta():
    theta = parse_scalar("theta")

    assert exp_plus_i_half(theta) == sp.exp(sp.I * theta / 2)


def test_numeric_float_angles_use_fast_float_coefficients():
    theta = sp.Float("0.78539816339744839", 17)

    assert isinstance(cos_half(theta), float)
    assert isinstance(sin_half(theta), float)
    assert abs(cos_half(theta) - 0.9238795325112867) < 1e-15
    assert abs(sin_half(theta) - 0.3826834323650898) < 1e-15


def test_scalar_simplify_skips_exact_simplify_for_float_numbers():
    value = scalar_simplify(sp.cos(sp.Float("0.78539816339744839", 17)))

    assert isinstance(value, float)


def test_high_precision_numeric_angles_become_rigorous_balls():
    theta = sp.Float(
        "0.78539816339744830961566084581987572104929234984377645524373614807695410157155225",
        80,
    )

    value = cos_half(theta)

    assert isinstance(value, acb)
    with ctx.workdps(80):
        expected = (to_acb(theta) / 2).cos()
        assert value.overlaps(expected)
        assert value.rad() < arb("1e-75")


def test_float_configuration_values_keep_their_decimal_spelling():
    with ctx.workdps(60):
        assert (to_arb(1e-4) - arb("1/10000")).contains(0)
