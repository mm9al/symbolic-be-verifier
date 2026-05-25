from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Dict, Mapping

import sympy as sp
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application,
    parse_expr,
    standard_transformations,
)


PAULI_ORDER = ("I", "X", "Y", "Z")
_IDENTITY_SYMBOL = sp.Symbol("Iop", commutative=False)
_X_SYMBOL = sp.Symbol("Xop", commutative=False)
_Y_SYMBOL = sp.Symbol("Yop", commutative=False)
_Z_SYMBOL = sp.Symbol("Zop", commutative=False)

_PAULI_PRODUCT = {
    ("I", "I"): (1, "I"),
    ("I", "X"): (1, "X"),
    ("I", "Y"): (1, "Y"),
    ("I", "Z"): (1, "Z"),
    ("X", "I"): (1, "X"),
    ("Y", "I"): (1, "Y"),
    ("Z", "I"): (1, "Z"),
    ("X", "X"): (1, "I"),
    ("Y", "Y"): (1, "I"),
    ("Z", "Z"): (1, "I"),
    ("X", "Y"): (sp.I, "Z"),
    ("Y", "Z"): (sp.I, "X"),
    ("Z", "X"): (sp.I, "Y"),
    ("Y", "X"): (-sp.I, "Z"),
    ("Z", "Y"): (-sp.I, "X"),
    ("X", "Z"): (-sp.I, "Y"),
}


def _clean_terms(terms: Mapping[str, sp.Expr]) -> Dict[str, sp.Expr]:
    cleaned: Dict[str, sp.Expr] = {}
    for pauli, coeff in terms.items():
        if pauli not in PAULI_ORDER:
            raise ValueError(f"Unsupported Pauli operator: {pauli!r}")
        simplified = sp.simplify(coeff)
        if simplified != 0:
            cleaned[pauli] = simplified
    return cleaned


@dataclass(frozen=True)
class OpExpr:
    """A symbolic operator over one system qubit in the Pauli basis."""

    terms: Mapping[str, sp.Expr]

    def __post_init__(self) -> None:
        object.__setattr__(self, "terms", _clean_terms(self.terms))

    @staticmethod
    def zero() -> "OpExpr":
        return OpExpr({})

    @staticmethod
    def identity() -> "OpExpr":
        return OpExpr({"I": sp.Integer(1)})

    @staticmethod
    def pauli(pauli: str) -> "OpExpr":
        if pauli not in PAULI_ORDER:
            raise ValueError(f"Unsupported Pauli operator: {pauli!r}")
        return OpExpr({pauli: sp.Integer(1)})

    def is_zero(self) -> bool:
        return not self.terms

    def simplify(self) -> "OpExpr":
        return OpExpr(self.terms)

    def scale(self, coeff: sp.Expr) -> "OpExpr":
        coeff = sp.sympify(coeff)
        return OpExpr({pauli: coeff * value for pauli, value in self.terms.items()})

    def __add__(self, other: "OpExpr") -> "OpExpr":
        terms = dict(self.terms)
        for pauli, coeff in other.terms.items():
            terms[pauli] = terms.get(pauli, sp.Integer(0)) + coeff
        return OpExpr(terms)

    def __sub__(self, other: "OpExpr") -> "OpExpr":
        return self + other.scale(-1)

    def __neg__(self) -> "OpExpr":
        return self.scale(-1)

    def __mul__(self, other: "OpExpr") -> "OpExpr":
        terms: Dict[str, sp.Expr] = {}
        for left_pauli, left_coeff in self.terms.items():
            for right_pauli, right_coeff in other.terms.items():
                phase, product = _PAULI_PRODUCT[(left_pauli, right_pauli)]
                terms[product] = terms.get(product, sp.Integer(0)) + left_coeff * right_coeff * phase
        return OpExpr(terms)

    def equals(self, other: "OpExpr") -> bool:
        diff = self - other
        return all(sp.simplify(coeff) == 0 for coeff in diff.terms.values())

    def __str__(self) -> str:
        if not self.terms:
            return "0"

        ordered = [(pauli, self.terms[pauli]) for pauli in PAULI_ORDER if pauli in self.terms]
        common = ordered[0][1]
        if len(ordered) > 1 and all(sp.simplify(coeff - common) == 0 for _, coeff in ordered):
            inside = " + ".join(pauli for pauli, _ in ordered)
            return _format_common_factor(inside, common)

        pieces = [_format_term(pauli, coeff) for pauli, coeff in ordered]
        text = pieces[0]
        for piece in pieces[1:]:
            if piece.startswith("-"):
                text += f" - {piece[1:]}"
            else:
                text += f" + {piece}"
        return text

    def __repr__(self) -> str:
        return f"OpExpr({self})"


def _format_scalar(coeff: sp.Expr) -> str:
    return sp.sstr(sp.simplify(coeff))


def _format_common_factor(inside: str, coeff: sp.Expr) -> str:
    coeff = sp.simplify(coeff)
    if coeff == 1:
        return inside
    if coeff == -1:
        return f"-({inside})"
    if coeff.is_Rational and coeff.p == 1:
        return f"({inside})/{coeff.q}"
    if coeff.is_Rational and coeff.p == -1:
        return f"-({inside})/{coeff.q}"
    return f"{_format_scalar(coeff)}*({inside})"


def _format_term(pauli: str, coeff: sp.Expr) -> str:
    coeff = sp.simplify(coeff)
    if coeff == 1:
        return pauli
    if coeff == -1:
        return f"-{pauli}"
    if coeff.is_Rational and coeff.p == 1:
        return f"{pauli}/{coeff.q}"
    if coeff.is_Rational and coeff.p == -1:
        return f"-{pauli}/{coeff.q}"
    return f"{_format_scalar(coeff)}*{pauli}"


def zero() -> OpExpr:
    return OpExpr.zero()


def identity() -> OpExpr:
    return OpExpr.identity()


def pauli(name: str) -> OpExpr:
    return OpExpr.pauli(name)


def parse_operator_expression(text: str) -> OpExpr:
    """Parse a small Pauli-basis expression such as ``(X + Z)/2``."""

    text = _preprocess_expression(text)
    local_dict = {
        "I": _IDENTITY_SYMBOL,
        "X": _X_SYMBOL,
        "Y": _Y_SYMBOL,
        "Z": _Z_SYMBOL,
        "i": sp.I,
        "pi": sp.pi,
        "sqrt": sp.sqrt,
    }
    transformations = standard_transformations + (implicit_multiplication_application,)
    parsed = parse_expr(text, local_dict=local_dict, transformations=transformations, evaluate=True)
    return _sympy_to_op(parsed)


def _preprocess_expression(text: str) -> str:
    stripped = text.strip()
    stripped = stripped.replace("^", "**")
    stripped = re.sub(r"\bi(?=[IXYZ(])", "i*", stripped)
    return stripped


def _sympy_to_op(expr: sp.Expr) -> OpExpr:
    expr = sp.expand(expr)
    result = OpExpr.zero()
    for term in sp.Add.make_args(expr):
        result += _sympy_term_to_op(term)
    return result


def _sympy_term_to_op(term: sp.Expr) -> OpExpr:
    symbols = {
        _IDENTITY_SYMBOL: "I",
        _X_SYMBOL: "X",
        _Y_SYMBOL: "Y",
        _Z_SYMBOL: "Z",
    }

    coeff = sp.Integer(1)
    op = OpExpr.identity()
    for factor in sp.Mul.make_args(term):
        if factor in symbols:
            op = op * OpExpr.pauli(symbols[factor])
        elif isinstance(factor, sp.Pow) and factor.base in symbols and factor.exp.is_integer:
            exponent = int(factor.exp)
            if exponent < 0:
                raise ValueError(f"Negative Pauli powers are not supported: {factor}")
            for _ in range(exponent):
                op = op * OpExpr.pauli(symbols[factor.base])
        else:
            coeff *= factor
    return op.scale(coeff)
