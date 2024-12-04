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

        I_elements = [
                [[kronecker_delta(i, k) * sum([kronecker_delta(i, l) * (self.elements[l-1][l-1][l-1])**(-1) for l in range(1, self.n + 1)]) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]


        J_elements = [
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
                [[kronecker_delta(i, j) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ]

        return (Cubriz(I_elements), Cubriz(J_elements))

    def inverse_IAJ_jk(self):
        assert self.m == self.n == self.o

        I_elements = [
                [[kronecker_delta(j, k) * sum([kronecker_delta(j, l) * (self.elements[l-1][l-1][l-1])**(-1) for l in range(1, self.n + 1)]) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
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
                [[kronecker_delta(i, j) * sum([kronecker_delta(i, l) * (self.elements[l-1][l-1][l-1])**(-1) for l in range(1, self.n + 1)]) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
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

A2 = Cubriz([
    [[6.0, 8.0],
     [4.0, 2.0]],

    [[7.0, 6.0],
     [5.0, 4.0]]
])

A = Cubriz([
    [[6.0, 8.0, 4.0],
     [2.0, 7.0, 7.0],
     [2.0, 7.0, 1.0]],

    [[7.0, 6.0, 5.0],
     [4.0, 3.0, 2.0],
     [1.0, 10.0, 9.0]],

    [[1.0, 2.0, 3.0],
     [6.0, 5.0, 4.0],
     [7.0, 80.0, 9.0]]
])

def print_expansion(A):
    for i in range(A.m):
        for j in range(A.n):
            for k in range(A.o):
                product = f""
                for l in range(A.n):
                    product += f" + A_{i}{l}{k} I_{l}{j}{k} J_{i}{j}{l}"

                print(f"B_{i}{j}{k} = delta_{i}{k} = {int(i == k)} =", product)


(I1, J1) = A.inverse_AIJ_ij()
(I2, J2) = A.inverse_AIJ_jk()
(I3, J3) = A.inverse_AIJ_ik()

(I4, J4) = A.inverse_IAJ_ij()
(I5, J5) = A.inverse_IAJ_jk()
(I6, J6) = A.inverse_IAJ_ik()

(I7, J7) = A.inverse_IJA_ij()
(I8, J8) = A.inverse_IJA_jk()
(I9, J9) = A.inverse_IJA_ik()

print("Pares inversa sobre ij:")
print("AIJ:")
A.times(I1, J1).print()
print("IAJ:")
I4.times(A, J4).print()
print("IJA:")
I7.times(J7, A).print()

print("Pares inversa sobre jk:")
print("AIJ:")
A.times(I2, J2).print()
print("IAJ:")
I5.times(A, J5).print()
print("IJA:")
I8.times(J8, A).print()

print("Pares inversa sobre ik:")
print("AIJ:")
A.times(I3, J3).print()
print("IAJ:")
I6.times(A, J6).print()
print("IJA:")
I9.times(J9, A).print()
