from typing import Callable
from random import uniform

# M_oxmxn
# Slice  Row   Column
# A[k]   [i]   [j]
class Cubrix:
    o: int
    m: int
    n: int
    elements: list[list[list[float]]]

    def __init__(self, elements: list[list[list[float]]]):
        self.elements = elements
        self.o = len(elements)
        self.m = len(elements[0])
        self.n = len(elements[0][0])

    def from_function(n: int, func: Callable[[int, int, int], float]):
        return Cubrix([
            [[func(i, j, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
        ])

    def __add__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            return Cubrix([
                [[self.elements[k-1][i-1][j-1] + other for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
            ])
        elif isinstance(other, Cubrix):
            return Cubrix([
                [[self.elements[k-1][i-1][j-1] + other.elements[k-1][i-1][j-1] for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
            ])
        else:
            pass

    def __mul__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            return Cubrix([
                [[other * self.elements[k-1][i-1][j-1] for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
            ]) 
        else:
            pass

    def __pow__(self, exponent):
        return Cubrix([
                [[(self.elements[k-1][i-1][j-1])**(exponent) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ])

        return Cubrix(elements)

    def print(self, decimal_positions):
        s = ""
        for i in range(1, self.m + 1):
            s += "| "
            for k in range(1, self.o + 1):
                for j in range(1, self.n + 1):
                    s += f"{self.elements[k-1][i-1][j-1]:.{decimal_positions}f}"
                    s += " "
                s += "|"
            s += "\n"

        print(s)

    def inverse_pair(self):
        return (self**(-1), D(self.n))

    __radd__ = __add__
    __rmul__ = __mul__

    def _times_ijk(self, B, C, i, j, k):
        c = 0
        for l in range(1, self.n + 1):
            c += self.elements[k-1][i-1][l-1] * B.elements[k-1][l-1][j-1] * C.elements[l-1][i-1][j-1]
        return c

    def times(self, B, C):
        assert self.n == B.m
        assert self.n == C.o

        result_m = min(self.m, C.m)
        result_n = min(B.n, C.n)
        result_o = min(self.o, B.o)

        elements = [
            [[self._times_ijk(B, C, i, j, k) for j in range(1, result_n + 1)] for i in range(1, result_m + 1)]
            for k in range(1, result_o + 1)
        ]

        return Cubrix(elements)

def d(a, b):
    return float(a == b)
    
def D_ij(n):
    return Cubrix.from_function(n, lambda i, j, k: d(i, j))

def D_jk(n):
    return Cubrix.from_function(n, lambda i, j, k: d(j, k))

def D_ik(n):
    return Cubrix.from_function(n, lambda i, j, k: d(i, k))

def D(n):
    return Cubrix.from_function(n, lambda i, j, k: d(i, j) * d(j, k))

def I(n):
    return Cubrix.from_function(n, lambda i, j, k: 1.0)

def random(n, lower, upper):
    return Cubrix.from_function(n, lambda i, j, k: uniform(lower, upper))
