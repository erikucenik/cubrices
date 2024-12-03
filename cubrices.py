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

    def inverse_AIJ_ij(self):
        assert self.m == self.n == self.o

        I_elements = [
                [[kronecker_delta(i, j) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        J_elements = [
                [[kronecker_delta(i, j) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

    def inverse_AIJ_jk(self):
        assert self.m == self.n == self.o

        J_elements = [
                [[kronecker_delta(j, k) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        I_elements = [
                [[kronecker_delta(j, k) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

    def inverse_AIJ_ik(self):
        assert self.m == self.n == self.o

        J_elements = [
                [[kronecker_delta(i, k) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        I_elements = [
                [[kronecker_delta(i, k) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))
    

    def inverse_IAJ_ij(self):
        assert self.m == self.n == self.o

        I_elements = [
                [[kronecker_delta(i, j) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        J_elements = [
                [[kronecker_delta(i, j) * kronecker_delta(i, k) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

    def inverse_IAJ_jk(self):
        assert self.m == self.n == self.o

        I_elements = [
                [[kronecker_delta(j, k) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        J_elements = [
                [[kronecker_delta(j, k) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

    def inverse_IAJ_ik(self):
        assert self.m == self.n == self.o

        J_elements = [
                [[kronecker_delta(i, k) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        I_elements = [
                [[kronecker_delta(i, k) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

    def inverse_IJA_ij(self):
        assert self.m == self.n == self.o

        I_elements = [
                [[kronecker_delta(i, j) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        J_elements = [
                [[kronecker_delta(i, j) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

    def inverse_IJA_jk(self):
        assert self.m == self.n == self.o

        I_elements = [
                [[kronecker_delta(j, k) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        J_elements = [
                [[kronecker_delta(j, k) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

    def inverse_IJA_ik(self):
        assert self.m == self.n == self.o

        J_elements = [
                [[kronecker_delta(i, k) * (self.elements[k-1][i-1][j-1])**(-1) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        I_elements = [
                [[kronecker_delta(i, k) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

def kronecker_delta(a, b):
    return float(a == b)

def d_ij_cubrix(n):
    return Cubriz([
            [[kronecker_delta(i, j) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

A = Cubriz([
    [[6.0, 8.0],
     [4.0, 2.0]],

    [[7.0, 6.0],
     [5.0, 4.0]]
])

print("A:")
A.print()

(I_ij1, J_ij1) = A.inverse_AIJ_ij()
(I_jk1, J_jk1) = A.inverse_AIJ_jk()
(I_ik1, J_ik1) = A.inverse_AIJ_ik()

print("AIJ sobre ij")
A.times(I_ij1, J_ij1).print()

print("AIJ sobre jk")
A.times(I_jk1, J_jk1).print()

print("AIJ sobre ik (NO FUNCIONA)")
A.times(I_ik1, J_ik1).print()

(I_ij2, J_ij2) = A.inverse_IAJ_ij()
(I_jk2, J_jk2) = A.inverse_IAJ_jk()
(I_ik2, J_ik2) = A.inverse_IAJ_ik()

print("IAJ sobre ij")
I_ij2.times(A, J_ij2).print()
print("IAJ sobre jk (NO FUNCIONA)")
I_jk2.times(A, J_jk2).print()
print("IAJ sobre ik")
I_ik2.times(A, J_ik2).print()

(I_ij3, J_ij3) = A.inverse_IJA_ij()
(I_jk3, J_jk3) = A.inverse_IJA_jk()
(I_ik3, J_ik3) = A.inverse_IJA_ik()

print("IJA sobre ij (NO FUNCIONA)")
I_ij3.times(J_ij3, A).print()
print("IJA sobre jk")
I_jk3.times(J_jk3, A).print()
print("IJA sobre ik")
I_ik3.times(J_ik3, A).print()
