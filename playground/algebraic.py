import math
from collections import deque
from fractions import Fraction
import sys


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


def uspensky_positive(p: Polynomial, max_depth: int) -> list[Interval]:
    """Implementation of Uspensky's algorithm

    TODO: add citations
    """
    intervals: list[Interval] = []
    queue = deque(
        [
            (
                p,
                Fraction(0),
                cauchy_upper_bound(p),
                0,
            )
        ]
    )

    while len(queue) > 0:
        print("\n== Iteration ==")
        p, lower, upper, depth = queue.pop()
        print(f"= Depth {depth} =")
        print(f"Checking for root between {float(lower)} and {float(upper)}")
        print(p.coefs)
        print(intervals)

        # Test p
        v = var(p.coefs)
        if v == 0:
            # No roots contained in the [lower, upper] interval
            print("No roots contained in interval")
            continue
        if v == 1:
            # Exactly 1 root contained in the [lower, upper] interval
            print("1 root contained in interval")
            if depth == max_depth:
                intervals.append((lower, upper))
                print("terminating branch")
                continue
            print("searching further")

        if depth == max_depth:
            print("D:")
            continue

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

        # Sub 1 / (x + 1)
        coefs_b = []
        for i in range(len(p.coefs)):
            coef = Fraction(0)
            for j in range(i, len(p.coefs)):
                coef += Fraction(math.comb(j, i)) * p.coefs[len(p.coefs) - 1 - j]
            coefs_b.append(coef)
        b = Polynomial(coefs_b)

        lower_b = 1 / (upper + 1)
        upper_b = 1 / (lower + 1)
        print("bounds for b:", lower_b, upper_b)

        queue.append((b, lower_b, upper_b, depth + 1))
        queue.append((a, lower + 1, upper + 1, depth + 1))

    return intervals


p = Polynomial([Fraction(1), Fraction(0), Fraction(-3), Fraction(1)])


for interval in uspensky_positive(p, int(sys.argv[1])):
    print()
    avg = interval[0] / 2 + interval[1] / 2
    print(f"avg: {avg.numerator / avg.denominator}")
