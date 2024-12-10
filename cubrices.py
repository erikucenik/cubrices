from typing import Callable
from random import uniform

class Cubrix:
    """
    A class used to represent Cubrices: mathematical objects described
    in Erik's paper (see https://github.com/erikucenik/cubrices-paper).

    m, n, o : int
        The dimensions of the Cubrix. They must lie in [1, k].
        m is the number of rows.
        n is the number of columns.
        o is the number of slices.

        They are automatically determined in __init__.

    elements : list[list[list[float]]]
        The elements of the Cubrix.

        Be aware of the indexing. By definition, A.elements[a][b][c]
        will give you the element in the (a + 1)th slice, the
        (b + 1)th row and the (c + 1)th column. If you wish to
        get A_{ijk}, you must write A.elements[k-1][i-1][j-1].

        This design choice lies in practicality.
    """

    o: int
    m: int
    n: int
    elements: list[list[list[float]]]

    def __init__(self, elements: list[list[list[float]]]) -> None:
        self.elements = elements
        self.o = len(elements)
        self.m = len(elements[0])
        self.n = len(elements[0][0])

    def from_function(n: int, func: Callable[[int, int, int], float]) -> "Cubrix":
        """
        Returns a cubic Cubrix of size n whose elements are created
        using the return value obtained by calling func(i, j, k)
        on each i, j and k.
        """

        return Cubrix([
            [[func(i, j, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
        ])

    def __add__(self, other) -> "Cubrix":
        """
        Addition operation overloading. Let self be a Cubrix and Z be the output Cubrix.

        If other : float | int
            Z_{ijk} = self_{ijk} + other
        If other : Cubrix
            Z_{ijk} = self_{ijk} + other_{ijk}
            They must have equal dimensions. 
        """

        if isinstance(other, float) or isinstance(other, int):
            return Cubrix([
                [[self.elements[k-1][i-1][j-1] + other for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
            ])

        elif isinstance(other, Cubrix):
            assert self.m == other.m
            assert self.n == other.n
            assert self.o == other.o

            return Cubrix([
                [[self.elements[k-1][i-1][j-1] + other.elements[k-1][i-1][j-1] for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
            ])

        else:
            raise TypeError("Addition is not defined for types other than floats or other Cubrices.")

    def __mul__(self, scalar: float | int) -> "Cubrix":
        """
        Multiplication operation overloading. Scales each element of
        the Cubrix by the value scalar.

        Note: for Cubrix multiplication, see the times() method.
        """
        
        return Cubrix([
            [[scalar * self.elements[k-1][i-1][j-1] for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
            for k in range(1, self.o + 1)
        ]) 

    def __pow__(self, exponent: float | int) -> "Cubrix":
        """
        Power operation overloading. Raises every element of the Cubrix
        the the value exponent.
        """

        return Cubrix([
                [[(self.elements[k-1][i-1][j-1])**(exponent) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ])

    def to_str(self, decimal_positions: int) -> str:
        """
        Returns the two-dimensional representation of the Cubrix with
        its elements rounded up to as many positions as specified by
        the decimal_positions parameter.
        """

        s = ""
        for i in range(1, self.m + 1):
            s += "| "
            for k in range(1, self.o + 1):
                for j in range(1, self.n + 1):
                    s += f"{self.elements[k-1][i-1][j-1]:.{decimal_positions}f}"
                    s += " "
                s += "|"
            s += "\n"

        return s

    def print(self, decimal_positions: int) -> None:
        """
        Prints the Cubrix string representation rounded up to
        as many positions as specified by decimal_positions.
        """

        print(self.to_str(decimal_positions))

    def inverse_pair(self) -> tuple["Cubrix", "Cubrix"]:
        """
        Returns the inverse_pair (A^{-1}, \Delta) of self.
        The user will decide upon the order in which to do
        the product to get the corresponding Kronecker Cubrix.
        """

        return (self**(-1), D(self.n))

    __radd__ = __add__
    __rmul__ = __mul__

    def _times_ijk(self, B: "Cubrix", C: "Cubrix", i: int, j: int, k: int) -> float:
        """
        Underlying method for the times() method. Returns the (i,j,k)th
        value of the cubrix that results from multiplying self, B
        and C as defined by the paper.
        """

        c = 0.0
        for l in range(1, self.n + 1):
            c += self.elements[k-1][i-1][l-1] * B.elements[k-1][l-1][j-1] * C.elements[l-1][i-1][j-1]
        return c

    def times(self, B: "Cubrix", C: "Cubrix") -> "Cubrix":
        """
        Returns the product of self, B and C. The same constraints are
        applied as in the paper.
        """

        assert self.n == B.m
        assert self.n == C.o

        result_m = min(self.m, C.m)
        result_n = min(B.n, C.n)
        result_o = min(self.o, B.o)

        return Cubrix([
            [[self._times_ijk(B, C, i, j, k) for j in range(1, result_n + 1)] for i in range(1, result_m + 1)]
            for k in range(1, result_o + 1)
        ])


def d(a: float | int, b: float | int) -> float:
    """Kronecker delta."""
    return float(a == b)
    
def D_ij(n: int) -> "Cubrix":
    """Generates a cubic Kronecker Cubrix of size n over i, j."""
    return Cubrix.from_function(n, lambda i, j, k: d(i, j))

def D_jk(n: int) -> "Cubrix":
    """Generates a cubic Kronecker Cubrix of size n over j, k."""
    return Cubrix.from_function(n, lambda i, j, k: d(j, k))

def D_ik(n: int) -> "Cubrix":
    """Generates a cubic Kronecker Cubrix of size n over i, k."""
    return Cubrix.from_function(n, lambda i, j, k: d(i, k))

def D(n: int) -> "Cubrix":
    """Generates a cubic Kronecker Cubrix of size n over i, j, k."""
    return Cubrix.from_function(n, lambda i, j, k: d(i, j) * d(j, k))

def I(n: int) -> "Cubrix":
    """Generates a cubic Cubrix of size n composed by just 1.0s."""
    return Cubrix.from_function(n, lambda i, j, k: 1.0)

def random(n: int, lower: float, upper: float) -> "Cubrix":
    """Generates a cubic Cubrix of size n composed of random floats
    ranging from lower to upper."""
    return Cubrix.from_function(n, lambda i, j, k: uniform(lower, upper))
