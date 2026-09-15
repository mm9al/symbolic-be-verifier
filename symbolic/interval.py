from __future__ import annotations

"""Rigorous real/complex ball helpers used at the numerical boundary.

The symbolic engine continues to use SymPy for exact expressions.  As soon as
an inexact numerical value is needed for verification, it is converted to an
Arb/Acb ball and all subsequent numerical operations stay in ball arithmetic.
"""

from typing import Any

from flint import acb, arb, ctx
import sympy as sp


Ball = arb | acb


def is_ball(value: Any) -> bool:
    return isinstance(value, (arb, acb))


def is_exact_zero(value: Any) -> bool:
    if is_ball(value):
        return value.is_zero()
    return value == 0


def is_exact_one(value: Any) -> bool:
    return value == 1


def to_acb(value: Any) -> acb:
    """Evaluate a numeric SymPy expression as a rigorous complex ball."""

    if isinstance(value, acb):
        return value
    if isinstance(value, arb):
        return acb(value)
    if isinstance(value, bool):
        return acb(int(value))
    if isinstance(value, int):
        return acb(value)
    if isinstance(value, float):
        # Public float arguments are configuration values originating from
        # decimal CLI/CSV text (epsilon, tau, tolerances).  Reusing Python's
        # shortest round-trip decimal spelling preserves that intended paper
        # parameter instead of silently changing it to the binary64 ratio.
        return acb(arb(repr(value)))
    if isinstance(value, complex):
        return acb(to_arb(value.real), to_arb(value.imag))

    expr = sp.sympify(value)
    if expr == sp.I:
        return acb(0, 1)
    if isinstance(expr, sp.Integer):
        return acb(int(expr))
    if isinstance(expr, sp.Rational):
        return acb(arb(f"{expr.p}/{expr.q}"))
    if isinstance(expr, sp.Float):
        # A SymPy Float is an exact binary value.  Passing its mantissa and
        # exponent to Arb avoids an untracked decimal round trip.
        sign, mantissa, exponent, _ = expr._mpf_
        signed_mantissa = -mantissa if sign else mantissa
        return acb(arb((signed_mantissa, exponent)))
    if expr == sp.pi:
        return acb(arb.pi())
    if expr.is_Add:
        result = acb(0)
        for term in expr.args:
            result += to_acb(term)
        return result
    if expr.is_Mul:
        result = acb(1)
        for factor in expr.args:
            result *= to_acb(factor)
        return result
    if expr.is_Pow:
        base, exponent = expr.args
        if isinstance(exponent, sp.Integer):
            return to_acb(base) ** int(exponent)
        return to_acb(base) ** to_acb(exponent)
    if expr.func == sp.sin:
        return to_acb(expr.args[0]).sin()
    if expr.func == sp.cos:
        return to_acb(expr.args[0]).cos()
    if expr.func == sp.exp:
        return to_acb(expr.args[0]).exp()
    if expr.func == sp.conjugate:
        return to_acb(expr.args[0]).conjugate()
    if expr.func == sp.re:
        return acb(to_acb(expr.args[0]).real)
    if expr.func == sp.im:
        return acb(to_acb(expr.args[0]).imag)
    raise ValueError(f"Cannot rigorously evaluate non-numeric expression as Acb: {expr!r}")


def to_arb(value: Any) -> arb:
    ball = to_acb(value)
    if not ball.imag.is_zero():
        raise ValueError(f"Expected a real value, got {value!r}")
    return ball.real


def add(left: Any, right: Any) -> Any:
    if is_exact_zero(left):
        return right
    if is_exact_zero(right):
        return left
    if is_ball(left) or is_ball(right):
        return to_acb(left) + to_acb(right)
    return left + right


def subtract(left: Any, right: Any) -> Any:
    if is_exact_zero(right):
        return left
    if is_exact_zero(left):
        return -right
    if is_ball(left) or is_ball(right):
        return to_acb(left) - to_acb(right)
    return left - right


def multiply(left: Any, right: Any) -> Any:
    if is_exact_zero(left) or is_exact_zero(right):
        return acb(0) if is_ball(left) or is_ball(right) else 0
    if is_exact_one(left):
        return right
    if is_exact_one(right):
        return left
    if is_ball(left) or is_ball(right):
        return to_acb(left) * to_acb(right)
    return left * right


def divide(left: Any, right: Any) -> Any:
    if is_exact_one(right):
        return left
    if is_ball(left) or is_ball(right):
        return to_acb(left) / to_acb(right)
    return left / right


def format_ball(value: Ball, digits: int = 16) -> str:
    return value.str(digits).replace("j", "i")


__all__ = [
    "Ball",
    "acb",
    "add",
    "arb",
    "ctx",
    "divide",
    "format_ball",
    "is_ball",
    "is_exact_one",
    "is_exact_zero",
    "multiply",
    "subtract",
    "to_acb",
    "to_arb",
]
