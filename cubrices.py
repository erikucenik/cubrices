# M_oxmxn
#  Slice Filas Columnas
# m[k]   [i]   [j]
class Cubriz:
    o: int
    m: int
    n: int
    elements: list[list[list[float]]]

    def __init__(self, elements: list[list[list[float]]]):
        self.elements = elements
        self.o = len(elements)
        self.m = len(elements[0])
        self.n = len(elements[0][0])

    def print(self):
        for i in range(1, self.m + 1):
            print('| ', end='')
            for k in range(1, self.o + 1):
                for j in range(1, self.n + 1):
                    print(self.elements[k-1][i-1][j-1], end=' ')
                print('|', end='')
            print('\n', end='')

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
        return Cubriz(elements)

    def inverse_ij(self):
        assert self.m == self.n == self.o

        elements = [
                [[1.0 / (self.elements[k-1][i-1][j-1]) if i == j else 0 for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(elements), d_ij_cubrix(self.n))

    def inverse_ik(self):
        assert self.m == self.n == self.o

        elements = [
                [[1.0 / (self.elements[k-1][i-1][j-1]) if i == k else 0 for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(elements), d_ik_cubrix(self.n))

    def inverse_jk(self):
        assert self.m == self.n == self.o

        elements = [
                [[1.0 / (self.elements[k-1][i-1][j-1]) if j == k else 0 for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(elements), d_jk_cubrix(self.n))

def kronecker_delta(a, b):
    return float(a == b)

def d_ij_cubrix(n):
    return Cubriz([
            [[kronecker_delta(i, j) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

def d_jk_cubrix(n):
    return Cubriz([
            [[kronecker_delta(j, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

def d_ik_cubrix(n):
    return Cubriz([
            [[kronecker_delta(i, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

AIX = d_ij_cubrix(2)
IAX = d_ij_cubrix(2)
IXA = d_jk_cubrix(2)

AXJ = d_jk_cubrix(2)
XAJ = d_ik_cubrix(2)
XJA = d_ik_cubrix(2)

A = Cubriz([
    [[6.0, 8.0],
     [4.0, 2.0]],

    [[7.0, 6.0],
     [5.0, 4.0]]
])

B = Cubriz([
    [[20.0, 40.0],
     [80.0, 160.0]],

    [[320.0, 640.0],
     [640.0, 640.0]]
])

print("A:")
A.print()

(A_ij, d_ij) = A.inverse_ij()
(A_jk, d_jk) = A.inverse_jk()
(A_ik, d_ik) = A.inverse_ik()

print("IJ")

print("A * Aij * dij")
A.times(A_ij, d_ij).print()
print("A * dij * Aij")
A.times(d_ij, A_ij).print()

print("Aij * A * dij")
A_ij.times(A, d_ij).print()
print("Aij * dij * A")
A_ij.times(d_ij, A).print()

print("dij * Aij * A")
d_ij.times(A_ij, A).print()
print("dij * A * Aij")
d_ij.times(A, A_ij).print()

print("IK")
print("A * Aik * dik")
A.times(A_ik, d_ik).print()
print("A * dik * Aik")
A.times(d_ik, A_ik).print()
print("Aik * A * dik")
A_ik.times(A, d_ik).print()
print("Aik * dik * A")
A_ik.times(d_ik, A).print()
print("dik * A * Aik")
d_ik.times(A, A_ik).print()
print("dik * Aik * A")
d_ik.times(A_ik, A).print()

print("------")
print("JK")

print("A * Ajk * djk")
A.times(A_jk, d_jk).print()

print("***************** A * Aij * dij")
A.times(A_ij, d_ij).print()

print("A * djk * Ajk")
A.times(d_jk, A_jk).print()

print("Ajk * A * djk")
A_jk.times(A, d_jk).print()
print("Ajk * djk * A")
A_jk.times(d_jk, A).print()

print("djk * A * Ajk")
d_jk.times(A, A_jk).print()

print("djk * Ajk * A")
d_jk.times(A_jk, A).print()

print("------")

print("A * Aij * B")
A.times(A_ij, B).print()
