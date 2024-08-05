import math
import sys
from collections import deque
from fractions import Fraction


class Polynomial:
    coefs: list[Fraction]

    def __init__(self, coefs: list[Fraction]):
        self.coefs = coefs

    def eval(self, x: Fraction) -> Fraction:
        total = Fraction(0)
        for i, c in enumerate(self.coefs):
            total += c * x**i
        return total

    @staticmethod
    def fraction(x: Fraction) -> "Polynomial":
        """Returns minimal polynomial of x"""
        return Polynomial([-x, Fraction(1)])

    @staticmethod
    def sqrt(x: Fraction) -> "Polynomial":
        """Returns minimal polynomial of sqrt(x)"""
        return Polynomial([-x, Fraction(0), Fraction(1)])


class Algebraic:
    min_poly: Polynomial
    lower: Fraction
    upper: Fraction


def sign(x: Fraction) -> int:
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def var(b: list[Fraction]):
    b_prime = [x for x in b if x != 0]
    count = 0
    for i in range(len(b_prime) - 1):
        if sign(b_prime[i]) != sign(b_prime[i + 1]):
            count += 1
    return count


Interval = tuple[Fraction, Fraction]


def cauchy_upper_bound(p: Polynomial) -> Fraction:
    return 1 + max([abs(c) for c in p.coefs[:-1]]) / abs(p.coefs[-1])


def find_root_intervals_recursive(
    p: Polynomial,
    lower: Fraction,
    upper: Fraction,
    level: int,
) -> list[Interval]:
    """Implementation of Uspensky's algorithm"""
    v = var(p.coefs)
    if v == 0:
        # No roots contained in the [lower, upper] interval
        return []
    if v == 1:
        # Exactly 1 root contained in the [lower, upper] interval
        if level == 1:
            return [(lower, upper)]

    if level == 1:
        raise RuntimeError("More levels required")

    intervals: list[Interval] = []

    # Otherwise we transform the polynomial and try again
    # Sub x + 1
    coefs_a = []
    for i in range(len(p.coefs)):
        coef = Fraction(0)
        for j in range(i, len(p.coefs)):
            coef += Fraction(math.comb(j, i)) * p.coefs[j]
        coefs_a.append(coef)
    a = Polynomial(coefs_a)
    if coefs_a[0] == 0:
        intervals.append((lower + 1, lower + 1))

    for transformed_interval in find_root_intervals_recursive(
        a, lower, upper, level - 1
    ):
        l, u = transformed_interval
        interval = l + 1, u + 1
        intervals.append(interval)

    # Sub 1 / (x + 1)
    coefs_b = []
    for i in range(len(p.coefs)):
        coef = Fraction(0)
        for j in range(i, len(p.coefs)):
            coef += Fraction(math.comb(j, i)) * p.coefs[len(p.coefs) - 1 - j]
        coefs_b.append(coef)
    b = Polynomial(coefs_b)

    for transformed_interval in find_root_intervals_recursive(
        b, lower, upper, level - 1
    ):
        l, u = transformed_interval
        interval = 1 / (u + 1), 1 / (l + 1)
        intervals.append(interval)

    return intervals


def find_root_intervals(p: Polynomial, level: int) -> list[Interval]:
    intervals = []

    # Copy so we can be destructive
    p = Polynomial(p.coefs.copy())

    # Eliminate trivial zeroes
    has_zero_root = False
    while p.coefs[0] == 0:
        has_zero_root = True
        p.coefs.pop(0)

    if has_zero_root:
        intervals.append((Fraction(0), Fraction(0)))

    upper = cauchy_upper_bound(p)
    intervals += find_root_intervals_recursive(p, Fraction(0), upper, level)

    for i in range(len(p.coefs)):
        if i % 2 == 1:
            p.coefs[i] *= -1

    for interval in find_root_intervals_recursive(p, Fraction(0), upper, level):
        intervals.append((-interval[1], -interval[0]))

    return intervals


p = Polynomial([Fraction(1), Fraction(0), Fraction(-3), Fraction(1)])
# p = Polynomial([Fraction(-2), Fraction(0), Fraction(1)])


level = int(sys.argv[1])
for interval in find_root_intervals(p, level):
    print()
    avg = interval[0] / 2 + interval[1] / 2
    print(f"avg: {avg.numerator / avg.denominator}")
    print(f"exact: {avg}")
