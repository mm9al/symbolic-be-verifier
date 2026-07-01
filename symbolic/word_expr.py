from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, Iterable, Mapping, Tuple

import sympy as sp

from .scalar import scalar_simplify


Atom = str
Word = Tuple[Atom, ...]

ATOMS = ("H", "Hd", "G", "Gd", "A", "Ad", "C", "Cd")
GARBAGE_ATOMS = frozenset(("G", "Gd", "A", "Ad", "C", "Cd"))


@dataclass(frozen=True)
class WordExpr:
    """A noncommutative polynomial over QSP/QSVT block atoms."""

    terms: Mapping[Word, sp.Expr]

    def __post_init__(self) -> None:
        cleaned: Dict[Word, sp.Expr] = {}
        for word, coeff in self.terms.items():
            key = tuple(word)
            for atom in key:
                if atom not in ATOMS:
                    raise ValueError(f"Unsupported QSP atom: {atom!r}")
            simplified = scalar_simplify(coeff)
            if simplified != 0:
                cleaned[key] = simplified
        object.__setattr__(self, "terms", cleaned)

    @staticmethod
    def zero() -> "WordExpr":
        return WordExpr({})

    @staticmethod
    def identity() -> "WordExpr":
        return WordExpr({(): sp.Integer(1)})

    @staticmethod
    def atom(name: Atom) -> "WordExpr":
        if name not in ATOMS:
            raise ValueError(f"Unsupported QSP atom: {name!r}")
        return WordExpr({(name,): sp.Integer(1)})

    def is_zero(self) -> bool:
        return not self.terms

    def simplify(self) -> "WordExpr":
        return WordExpr(self.terms)

    def scale(self, coeff: sp.Expr) -> "WordExpr":
        coeff = sp.sympify(coeff)
        return WordExpr({word: coeff * value for word, value in self.terms.items()})

    def with_num_qubits(self, num_qubits: int) -> "WordExpr":
        return self

    def has_any_atom(self, atoms: Iterable[Atom]) -> bool:
        atom_set = set(atoms)
        return any(atom in atom_set for word in self.terms for atom in word)

    def replace_atom(self, old: Atom, new: Atom) -> "WordExpr":
        return WordExpr({tuple(new if atom == old else atom for atom in word): coeff for word, coeff in self.terms.items()})

    def __add__(self, other: "WordExpr") -> "WordExpr":
        terms = dict(self.terms)
        for word, coeff in other.terms.items():
            terms[word] = terms.get(word, sp.Integer(0)) + coeff
        return WordExpr(terms)

    def __sub__(self, other: "WordExpr") -> "WordExpr":
        return self + other.scale(-1)

    def __neg__(self) -> "WordExpr":
        return self.scale(-1)

    def __mul__(self, other: "WordExpr") -> "WordExpr":
        terms: Dict[Word, sp.Expr] = {}
        for left_word, left_coeff in self.terms.items():
            for right_word, right_coeff in other.terms.items():
                key = left_word + right_word
                terms[key] = terms.get(key, sp.Integer(0)) + left_coeff * right_coeff
        return WordExpr(terms)

    def equals(self, other: "WordExpr") -> bool:
        diff = self - other
        return all(scalar_simplify(coeff) == 0 for coeff in diff.terms.values())

    def __str__(self) -> str:
        if not self.terms:
            return "0"
        pieces = [_format_term(word, self.terms[word]) for word in sorted(self.terms, key=_sort_key)]
        text = pieces[0]
        for piece in pieces[1:]:
            if piece.startswith("-"):
                text += f" - {piece[1:]}"
            else:
                text += f" + {piece}"
        return text

    def __repr__(self) -> str:
        return f"WordExpr({self})"


_REWRITE_RULES: Mapping[Tuple[Atom, Atom], WordExpr] = {
    ("Ad", "A"): WordExpr.identity() - WordExpr({("Hd", "H"): sp.Integer(1)}),
    ("G", "Gd"): WordExpr.identity() - WordExpr({("H", "Hd"): sp.Integer(1)}),
    ("Ad", "C"): WordExpr({("Hd", "G"): sp.Integer(-1)}),
    ("Cd", "A"): WordExpr({("Gd", "H"): sp.Integer(-1)}),
    ("G", "Cd"): WordExpr({("H", "Ad"): sp.Integer(-1)}),
    ("C", "Gd"): WordExpr({("A", "Hd"): sp.Integer(-1)}),
    ("Gd", "G"): WordExpr.identity() - WordExpr({("Cd", "C"): sp.Integer(1)}),
    ("A", "Ad"): WordExpr.identity() - WordExpr({("C", "Cd"): sp.Integer(1)}),
}


@dataclass(frozen=True)
class _WordRewriteRule:
    replacement: Tuple[Tuple[Word, sp.Expr], ...]
    priority: int
    garbage_delta: int


_WORD_REWRITE_RULES: Mapping[Tuple[Atom, Atom], _WordRewriteRule] = {
    ("Ad", "A"): _WordRewriteRule((((), sp.Integer(1)), (("Hd", "H"), sp.Integer(-1))), 0, -2),
    ("G", "Gd"): _WordRewriteRule((((), sp.Integer(1)), (("H", "Hd"), sp.Integer(-1))), 0, -2),
    ("Ad", "C"): _WordRewriteRule(((("Hd", "G"), sp.Integer(-1)),), 1, -1),
    ("Cd", "A"): _WordRewriteRule(((("Gd", "H"), sp.Integer(-1)),), 1, -1),
    ("G", "Cd"): _WordRewriteRule(((("H", "Ad"), sp.Integer(-1)),), 1, -1),
    ("C", "Gd"): _WordRewriteRule(((("A", "Hd"), sp.Integer(-1)),), 1, -1),
    ("Gd", "G"): _WordRewriteRule((((), sp.Integer(1)), (("Cd", "C"), sp.Integer(-1))), 2, 0),
    ("A", "Ad"): _WordRewriteRule((((), sp.Integer(1)), (("C", "Cd"), sp.Integer(-1))), 2, 0),
}


def atom(name: Atom) -> WordExpr:
    return WordExpr.atom(name)


def zero() -> WordExpr:
    return WordExpr.zero()


def identity() -> WordExpr:
    return WordExpr.identity()


def rewrite_once(expr: WordExpr) -> tuple[WordExpr, bool]:
    for word, coeff in sorted(expr.terms.items(), key=lambda item: _sort_key(item[0])):
        for index in range(len(word) - 1):
            pair = (word[index], word[index + 1])
            replacement = _REWRITE_RULES.get(pair)
            if replacement is None:
                continue

            prefix = WordExpr({word[:index]: coeff})
            suffix = WordExpr({word[index + 2 :]: sp.Integer(1)})
            rewritten_term = prefix * replacement * suffix
            remaining = dict(expr.terms)
            del remaining[word]
            return WordExpr(remaining) + rewritten_term, True
    return expr, False


def eliminate(expr: WordExpr) -> WordExpr:
    terms: Dict[Word, sp.Expr] = defaultdict(lambda: sp.Integer(0))
    for word, coeff in expr.terms.items():
        for normalized_word, normalized_coeff in _normalize_word_items(tuple(word)):
            terms[normalized_word] += coeff * normalized_coeff
    return WordExpr(terms)


def normalize_word(word: Word) -> WordExpr:
    return WordExpr(dict(_normalize_word_items(tuple(word))))


@lru_cache(maxsize=None)
def _normalize_word_items(word: Word) -> Tuple[Tuple[Word, sp.Expr], ...]:
    redex = _find_best_redex(word)
    if redex is None:
        return ((word, sp.Integer(1)),)

    index, rule = redex
    terms: Dict[Word, sp.Expr] = defaultdict(lambda: sp.Integer(0))
    for replacement_word, replacement_coeff in rule.replacement:
        new_word = word[:index] + replacement_word + word[index + 2 :]
        for normalized_word, normalized_coeff in _normalize_word_items(new_word):
            terms[normalized_word] += replacement_coeff * normalized_coeff

    return _normalize_word_result_items(terms)


def _find_best_redex(word: Word) -> tuple[int, _WordRewriteRule] | None:
    candidates = []
    for index in range(len(word) - 1):
        pair = (word[index], word[index + 1])
        rule = _WORD_REWRITE_RULES.get(pair)
        if rule is not None:
            candidates.append((rule.priority, rule.garbage_delta, index, rule))

    if not candidates:
        return None

    _, _, index, rule = min(candidates, key=lambda item: (item[0], item[1], item[2]))
    return index, rule


def _normalize_word_result_items(terms: Mapping[Word, sp.Expr]) -> Tuple[Tuple[Word, sp.Expr], ...]:
    return tuple(
        sorted(
            ((word, coeff) for word, coeff in terms.items() if coeff != 0),
            key=lambda item: _sort_key(item[0]),
        )
    )


def _sort_key(word: Word) -> Tuple[int, Tuple[int, ...]]:
    rank = {atom_name: index for index, atom_name in enumerate(ATOMS)}
    return (len(word), tuple(rank[atom_name] for atom_name in word))


def _format_term(word: Word, coeff: sp.Expr) -> str:
    coeff = scalar_simplify(coeff)
    body = "I" if not word else " ".join(word)
    if coeff == 1:
        return body
    if coeff == -1:
        return f"-{body}"
    if coeff.is_Rational and coeff.p == 1:
        return f"{body}/{coeff.q}"
    if coeff.is_Rational and coeff.p == -1:
        return f"-{body}/{coeff.q}"
    scalar = sp.sstr(coeff).replace("I", "i")
    if coeff.is_Add:
        scalar = f"({scalar})"
    return f"{scalar}*{body}"
