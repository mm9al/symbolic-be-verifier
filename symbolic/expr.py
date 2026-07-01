from __future__ import annotations

from dataclasses import dataclass
import re
import time
from typing import Dict, Iterable, Mapping, Optional, Tuple

import sympy as sp
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application,
    parse_expr,
    rationalize,
    standard_transformations,
)

from . import profile
from .scalar import scalar_simplify


PAULI_ORDER = ("I", "X", "Y", "Z")
PauliString = Tuple[str, ...]

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


def _clean_terms(terms: Mapping[PauliString, sp.Expr], num_qubits: int) -> Dict[PauliString, sp.Expr]:
    cleaned: Dict[PauliString, sp.Expr] = {}
    for pauli_string, coeff in terms.items():
        key = tuple(pauli_string)
        if len(key) != num_qubits:
            raise ValueError(f"Pauli string {key!r} has length {len(key)}, expected {num_qubits}")
        for pauli_op in key:
            if pauli_op not in PAULI_ORDER:
                raise ValueError(f"Unsupported Pauli operator: {pauli_op!r}")
        simplified = scalar_simplify(coeff)
        if simplified != 0:
            cleaned[key] = simplified
    return cleaned


@dataclass(frozen=True)
class OpExpr:
    """A symbolic operator over system qubits in the Pauli-string basis."""

    terms: Mapping[PauliString, sp.Expr]
    num_qubits: int = 1

    def __post_init__(self) -> None:
        if self.num_qubits < 1:
            raise ValueError("num_qubits must be at least 1")
        object.__setattr__(self, "terms", _clean_terms(self.terms, self.num_qubits))

    @staticmethod
    def zero(num_qubits: int = 1) -> "OpExpr":
        return OpExpr({}, num_qubits=num_qubits)

    @staticmethod
    def identity(num_qubits: int = 1) -> "OpExpr":
        return OpExpr({_identity_key(num_qubits): sp.Integer(1)}, num_qubits=num_qubits)

    @staticmethod
    def pauli(pauli_op: str, index: int = 0, num_qubits: int = 1) -> "OpExpr":
        if pauli_op not in PAULI_ORDER:
            raise ValueError(f"Unsupported Pauli operator: {pauli_op!r}")
        if index < 0 or index >= num_qubits:
            raise ValueError(f"Pauli index {index} is outside a {num_qubits}-qubit system")
        key = ["I"] * num_qubits
        key[index] = pauli_op
        return OpExpr({tuple(key): sp.Integer(1)}, num_qubits=num_qubits)

    @staticmethod
    def pauli_string(ops: Iterable[str]) -> "OpExpr":
        key = tuple(ops)
        if not key:
            raise ValueError("Pauli string must contain at least one operator")
        return OpExpr({key: sp.Integer(1)}, num_qubits=len(key))

    def is_zero(self) -> bool:
        return not self.terms

    def simplify(self) -> "OpExpr":
        return OpExpr(self.terms, num_qubits=self.num_qubits)

    def scale(self, coeff: sp.Expr) -> "OpExpr":
        coeff = sp.sympify(coeff)
        return OpExpr(
            {pauli_string: coeff * value for pauli_string, value in self.terms.items()},
            num_qubits=self.num_qubits,
        )

    def with_num_qubits(self, num_qubits: int) -> "OpExpr":
        if num_qubits == self.num_qubits:
            return self
        if num_qubits < self.num_qubits:
            raise ValueError(f"Cannot shrink {self.num_qubits}-qubit operator to {num_qubits} qubits")
        padded_terms = {
            pauli_string + ("I",) * (num_qubits - self.num_qubits): coeff
            for pauli_string, coeff in self.terms.items()
        }
        return OpExpr(padded_terms, num_qubits=num_qubits)

    def __add__(self, other: "OpExpr") -> "OpExpr":
        left, right = _align_exprs(self, other)
        started = time.perf_counter() if profile.current() is not None else None
        try:
            terms = dict(left.terms)
            for pauli_string, coeff in right.terms.items():
                terms[pauli_string] = terms.get(pauli_string, sp.Integer(0)) + coeff
        finally:
            if started is not None:
                profile.add_combine(time.perf_counter() - started)
        return OpExpr(terms, num_qubits=left.num_qubits)

    def __sub__(self, other: "OpExpr") -> "OpExpr":
        return self + other.scale(-1)

    def __neg__(self) -> "OpExpr":
        return self.scale(-1)

    def __mul__(self, other: "OpExpr") -> "OpExpr":
        left, right = _align_exprs(self, other)
        started = time.perf_counter() if profile.current() is not None else None
        try:
            terms: Dict[PauliString, sp.Expr] = {}
            for left_string, left_coeff in left.terms.items():
                for right_string, right_coeff in right.terms.items():
                    phase = sp.Integer(1)
                    product = []
                    for left_op, right_op in zip(left_string, right_string):
                        local_phase, local_product = _PAULI_PRODUCT[(left_op, right_op)]
                        phase *= local_phase
                        product.append(local_product)
                    key = tuple(product)
                    terms[key] = terms.get(key, sp.Integer(0)) + left_coeff * right_coeff * phase
        finally:
            if started is not None:
                profile.add_combine(time.perf_counter() - started)
        return OpExpr(terms, num_qubits=left.num_qubits)

    def equals(self, other: "OpExpr") -> bool:
        diff = self - other
        return all(scalar_simplify(coeff) == 0 for coeff in diff.terms.values())

    def __str__(self) -> str:
        if not self.terms:
            return "0"

        ordered = [(key, self.terms[key]) for key in sorted(self.terms, key=_sort_key)]
        common = ordered[0][1]
        if len(ordered) > 1 and all(scalar_simplify(coeff - common) == 0 for _, coeff in ordered):
            inside = " + ".join(_format_pauli_string(key, self.num_qubits) for key, _ in ordered)
            return _format_common_factor(inside, common)

        pieces = [_format_term(key, coeff, self.num_qubits) for key, coeff in ordered]
        text = pieces[0]
        for piece in pieces[1:]:
            if piece.startswith("-"):
                text += f" - {piece[1:]}"
            else:
                text += f" + {piece}"
        return text

    def __repr__(self) -> str:
        return f"OpExpr({self})"


def _identity_key(num_qubits: int) -> PauliString:
    return ("I",) * num_qubits


def _align_exprs(left: OpExpr, right: OpExpr) -> Tuple[OpExpr, OpExpr]:
    num_qubits = max(left.num_qubits, right.num_qubits)
    return left.with_num_qubits(num_qubits), right.with_num_qubits(num_qubits)


def _sort_key(pauli_string: PauliString) -> Tuple[int, ...]:
    rank = {op: index for index, op in enumerate(PAULI_ORDER)}
    return tuple(rank[op] for op in pauli_string)


def _format_pauli_string(pauli_string: PauliString, num_qubits: int) -> str:
    active = [(index, op) for index, op in enumerate(pauli_string) if op != "I"]
    if not active:
        return "I"
    if num_qubits == 1:
        return active[0][1]
    return " ".join(f"{op}_{index}" for index, op in active)


def _format_scalar(coeff: sp.Expr) -> str:
    return sp.sstr(scalar_simplify(coeff)).replace("I", "i")


def _format_common_factor(inside: str, coeff: sp.Expr) -> str:
    coeff = scalar_simplify(coeff)
    coeff_expr = sp.sympify(coeff)
    if coeff == 1:
        return inside
    if coeff == -1:
        return f"-({inside})"
    if coeff_expr.is_Rational and coeff_expr.p == 1:
        return f"({inside})/{coeff_expr.q}"
    if coeff_expr.is_Rational and coeff_expr.p == -1:
        return f"-({inside})/{coeff_expr.q}"
    return f"{_format_scalar(coeff)}*({inside})"


def _format_term(pauli_string: PauliString, coeff: sp.Expr, num_qubits: int) -> str:
    coeff = scalar_simplify(coeff)
    coeff_expr = sp.sympify(coeff)
    body = _format_pauli_string(pauli_string, num_qubits)
    if coeff == 1:
        return body
    if coeff == -1:
        return f"-{body}"
    if coeff_expr.is_Rational and coeff_expr.p == 1:
        return f"{body}/{coeff_expr.q}"
    if coeff_expr.is_Rational and coeff_expr.p == -1:
        return f"-{body}/{coeff_expr.q}"
    return f"{_format_factor_scalar(coeff)}*{body}"


def _format_factor_scalar(coeff: sp.Expr) -> str:
    text = _format_scalar(coeff)
    if sp.sympify(coeff).is_Add:
        return f"({text})"
    return text


def zero(num_qubits: int = 1) -> OpExpr:
    return OpExpr.zero(num_qubits)


def identity(num_qubits: int = 1) -> OpExpr:
    return OpExpr.identity(num_qubits)


def pauli(name: str, index: int = 0, num_qubits: int = 1) -> OpExpr:
    return OpExpr.pauli(name, index=index, num_qubits=num_qubits)


def parse_operator_expression(text: str, num_qubits: Optional[int] = None) -> OpExpr:
    """Parse expressions such as ``(X + Z)/2`` or ``(X0 X1 + Z0 Z1)/2``."""

    text = _preprocess_expression(text)
    symbols, inferred_num_qubits = _make_local_symbols(text)
    local_dict = {name: symbol for name, symbol in symbols.items()}
    local_dict.update({"i": sp.I, "pi": sp.pi, "sqrt": sp.sqrt, "sin": sp.sin, "cos": sp.cos, "exp": sp.exp})
    transformations = standard_transformations + (implicit_multiplication_application, rationalize)
    parsed = parse_expr(text, local_dict=local_dict, transformations=transformations, evaluate=True)

    target_num_qubits = num_qubits or inferred_num_qubits or 1
    return _sympy_to_op(parsed, symbols, target_num_qubits)


def _preprocess_expression(text: str) -> str:
    stripped = text.strip()
    stripped = stripped.replace("^", "**")
    stripped = re.sub(r"\b([IXYZ])_(\d+)\b", r"\1\2", stripped)
    stripped = re.sub(r"\bi(?=[IXYZ(])", "i*", stripped)
    return stripped


def _make_local_symbols(text: str) -> Tuple[Dict[str, sp.Symbol], Optional[int]]:
    names = set(re.findall(r"\b[IXYZ]\d*\b", text))
    if not names:
        names = {"I", "X", "Y", "Z"}

    inferred_indices = []
    symbols: Dict[str, sp.Symbol] = {}
    for name in names:
        symbols[name] = sp.Symbol(f"{name}op", commutative=False)
        if len(name) > 1:
            inferred_indices.append(int(name[1:]))

    for base in PAULI_ORDER:
        symbols.setdefault(base, sp.Symbol(f"{base}op", commutative=False))

    inferred_num_qubits = max(inferred_indices) + 1 if inferred_indices else None
    return symbols, inferred_num_qubits


def _sympy_to_op(expr: sp.Expr, symbols: Mapping[str, sp.Symbol], num_qubits: int) -> OpExpr:
    expr = sp.expand(expr)
    result = OpExpr.zero(num_qubits)
    for term in sp.Add.make_args(expr):
        result += _sympy_term_to_op(term, symbols, num_qubits)
    return result


def _sympy_term_to_op(term: sp.Expr, symbols: Mapping[str, sp.Symbol], num_qubits: int) -> OpExpr:
    symbol_to_pauli = {symbol: _symbol_name_to_pauli(name, num_qubits) for name, symbol in symbols.items()}

    coeff = sp.Integer(1)
    op = OpExpr.identity(num_qubits)
    for factor in sp.Mul.make_args(term):
        if factor in symbol_to_pauli:
            op = op * symbol_to_pauli[factor]
        elif isinstance(factor, sp.Pow) and factor.base in symbol_to_pauli and factor.exp.is_integer:
            exponent = int(factor.exp)
            if exponent < 0:
                raise ValueError(f"Negative Pauli powers are not supported: {factor}")
            for _ in range(exponent):
                op = op * symbol_to_pauli[factor.base]
        else:
            coeff *= factor
    return op.scale(scalar_simplify(coeff))


def _symbol_name_to_pauli(name: str, num_qubits: int) -> OpExpr:
    op = name[0]
    if len(name) == 1:
        return OpExpr.pauli(op, index=0, num_qubits=num_qubits)
    return OpExpr.pauli(op, index=int(name[1:]), num_qubits=num_qubits)
