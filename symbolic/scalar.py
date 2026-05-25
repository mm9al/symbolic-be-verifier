from __future__ import annotations

import sympy as sp
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application,
    parse_expr,
    rationalize,
    standard_transformations,
)


def parse_scalar(text: str) -> sp.Expr:
    """Parse a scalar expression used in QASM parameters or coefficients."""

    local_dict = {
        "I": sp.I,
        "i": sp.I,
        "pi": sp.pi,
        "sqrt": sp.sqrt,
        "sin": sp.sin,
        "cos": sp.cos,
        "exp": sp.exp,
    }
    transformations = standard_transformations + (implicit_multiplication_application, rationalize)
    return scalar_simplify(parse_expr(text, local_dict=local_dict, transformations=transformations, evaluate=True))


def scalar_simplify(expr: sp.Expr) -> sp.Expr:
    """Simplify one scalar coefficient without touching Pauli-string structure."""

    simplified = sp.sympify(expr)
    simplified = sp.simplify(simplified)
    simplified = sp.trigsimp(simplified)
    return sp.simplify(simplified)


def cos_half(theta: sp.Expr) -> sp.Expr:
    return scalar_simplify(sp.cos(theta / 2))


def sin_half(theta: sp.Expr) -> sp.Expr:
    return scalar_simplify(sp.sin(theta / 2))


def exp_minus_i_half(theta: sp.Expr) -> sp.Expr:
    return scalar_simplify(sp.exp(-sp.I * theta / 2))


def exp_plus_i_half(theta: sp.Expr) -> sp.Expr:
    return scalar_simplify(sp.exp(sp.I * theta / 2))
