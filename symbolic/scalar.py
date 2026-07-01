from __future__ import annotations

import cmath
from functools import lru_cache
import math
import time
from numbers import Number

import sympy as sp
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application,
    parse_expr,
    rationalize,
    standard_transformations,
)

from . import profile


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

    if profile.current() is None:
        return _scalar_simplify_impl(expr)

    started = time.perf_counter()
    try:
        return _scalar_simplify_impl(expr)
    finally:
        profile.add_simplify(time.perf_counter() - started)


def _scalar_simplify_impl(expr: sp.Expr) -> sp.Expr:
    if isinstance(expr, Number) and not isinstance(expr, sp.Basic):
        return expr

    sympified = sp.sympify(expr)
    if getattr(sympified, "is_Float", False):
        return float(sympified)
    if isinstance(sympified, sp.Basic) and sympified.has(sp.Float):
        if sympified.is_number:
            value = complex(sympified.evalf())
            return value.real if value.imag == 0 else value
        return sympified

    return _scalar_simplify_cached(sympified)


@lru_cache(maxsize=1_000_000)
def _scalar_simplify_cached(expr: sp.Expr) -> sp.Expr:
    simplified = sp.cancel(expr)
    simplified = sp.expand_mul(simplified)
    simplified = sp.simplify(simplified)
    simplified = _rewrite_unit_trig_squares(simplified)
    simplified = sp.collect(simplified, sp.I, evaluate=True)
    return simplified


def _rewrite_unit_trig_squares(expr: sp.Expr) -> sp.Expr:
    if not expr.is_Add:
        return expr

    terms = list(expr.args)
    groups: dict[tuple[sp.Expr, sp.Expr], dict[str, list[int]]] = {}
    for index, term in enumerate(terms):
        split = _split_trig_square_term(term)
        if split is None:
            continue
        kind, arg, coeff = split
        groups.setdefault((coeff, arg), {"sin": [], "cos": []})[kind].append(index)

    matched: set[int] = set()
    replacements: list[sp.Expr] = []
    for (coeff, _), matches in groups.items():
        pair_count = min(len(matches["sin"]), len(matches["cos"]))
        if pair_count == 0:
            continue
        matched.update(matches["sin"][:pair_count])
        matched.update(matches["cos"][:pair_count])
        replacements.extend(coeff for _ in range(pair_count))

    if not replacements:
        return expr

    kept = [term for index, term in enumerate(terms) if index not in matched]
    return sp.Add(*kept, *replacements)


def _split_trig_square_term(term: sp.Expr) -> tuple[str, sp.Expr, sp.Expr] | None:
    if term.is_Mul:
        trig_factor = None
        coeff_factors = []
        for factor in term.args:
            split = _split_bare_trig_square(factor)
            if split is None:
                coeff_factors.append(factor)
                continue
            if trig_factor is not None:
                return None
            trig_factor = split
        if trig_factor is None:
            return None
        kind, arg = trig_factor
        return kind, arg, sp.Mul(*coeff_factors)

    split = _split_bare_trig_square(term)
    if split is None:
        return None
    kind, arg = split
    return kind, arg, sp.Integer(1)


def _split_bare_trig_square(expr: sp.Expr) -> tuple[str, sp.Expr] | None:
    if not expr.is_Pow or expr.exp != 2:
        return None
    base = expr.base
    if base.func == sp.sin:
        return "sin", base.args[0]
    if base.func == sp.cos:
        return "cos", base.args[0]
    return None


def cos_half(theta: sp.Expr) -> sp.Expr:
    numeric = _numeric_float(theta)
    if numeric is not None:
        return math.cos(numeric / 2)
    return scalar_simplify(sp.cos(theta / 2))


def sin_half(theta: sp.Expr) -> sp.Expr:
    numeric = _numeric_float(theta)
    if numeric is not None:
        return math.sin(numeric / 2)
    return scalar_simplify(sp.sin(theta / 2))


def exp_minus_i_half(theta: sp.Expr) -> sp.Expr:
    numeric = _numeric_float(theta)
    if numeric is not None:
        return cmath.exp(-0.5j * numeric)
    return scalar_simplify(sp.exp(-sp.I * theta / 2))


def exp_plus_i_half(theta: sp.Expr) -> sp.Expr:
    numeric = _numeric_float(theta)
    if numeric is not None:
        return cmath.exp(0.5j * numeric)
    return scalar_simplify(sp.exp(sp.I * theta / 2))


def _numeric_float(theta: sp.Expr) -> float | None:
    if isinstance(theta, float):
        return theta
    if isinstance(theta, sp.Basic) and theta.has(sp.Float) and theta.is_number:
        return float(theta.evalf())
    return None
