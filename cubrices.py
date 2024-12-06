# M_oxmxn
#  Slice Filas Columnas
# A[k]   [i]   [j]
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

    def __add__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            return Cubriz([
                [[self.elements[k-1][i-1][j-1] + other for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
            ])
        elif isinstance(other, Cubriz):
            return Cubriz([
                [[self.elements[k-1][i-1][j-1] + other.elements[k-1][i-1][j-1] for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
            ])
        else:
            pass

    def __mul__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            return Cubriz([
                [[other * self.elements[k-1][i-1][j-1] for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
            ]) 
        else:
            pass

    def __pow__(self, exponent):
        return Cubriz([
                [[(self.elements[k-1][i-1][j-1])**(exponent) for j in range(1, self.n + 1)] for i in range(1, self.m + 1)]
                for k in range(1, self.o + 1)
        ])

        return Cubriz(elements)

    def inverse_pair(self):
        return (self**(-1), D(self.n))

    __radd__ = __add__
    __rmul__ = __mul__

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

def kronecker_delta(a, b):
    return float(a == b)

def d(a, b):
    return kronecker_delta(a, b)

def I(n):
    return Cubriz([
            [[1.0] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])
    
def D(n):
    return Cubriz([
        [[d(i, j) * d(j, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
        for k in range(1, n + 1)
    ])

def d_ij_cubrix(n):
    return Cubriz([
            [[kronecker_delta(i, j) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

def d_ik_cubrix(n):
    return Cubriz([
            [[kronecker_delta(i, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

def d_jk_cubrix(n):
    return Cubriz([
            [[kronecker_delta(j, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

A = Cubriz([
    [[6.0, 8.0],
     [4.0, 2.0]],

    [[7.0, 6.0],
     [5.0, 4.0]]
])

A3 = Cubriz([
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
def D(n):
    return Cubriz([
            [[d(i, j) * d(j, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

def I(n):
    return Cubriz([
            [[1 for j in range(1, n + 1)] for i in range(1, n + 1)]
            for k in range(1, n + 1)
    ])

def AI(A):
    return Cubriz([
        [[A.elements[k-1][k-1][k-1]**(-1) for j in range(1, A.n + 1)] for i in range(1, A.n + 1)]
        for k in range(1, A.n + 1)
    ])

#def AJ(n):
#    return Cubriz([
#        [[d(i, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
#        for k in range(1, n + 1)
#    ])

def AJ(n):
    return Cubriz([
        [[d(i, k) for j in range(1, n + 1)] for i in range(1, n + 1)]
        for k in range(1, n + 1)
    ])

"""
for i in range(1, A.n + 1):
    for j in range(1, A.n + 1):
        for k in range(1, A.n + 1):
            txt = f"B_{i}{j}{k} = d_{i}{k} = {d(i, k)} = "
            for l in range(1, A.n + 1):
                txt += f"+ A_{i}{l}{k} (d_{l}{k} A_{l}{l}{l}^-1) d_{i}{l} "

            print(txt)
for i in range(1, A.n + 1):
    for j in range(1, A.n + 1):
        for k in range(1, A.n + 1):
            txt = f"d_{i}{k} = {d(i, k)} = "
            for l in range(1, A.n + 1):
                #txt += f"+ A_{i}{l}{k} A_{l}{l}{l}^-1 (d_{l}{k} d_{i}{l})"
                txt += f"+ A_{i}{l}{k} A_{l}{l}{l}^-1 (d_{l}{k} d_{i}{l})"

            print(txt)

A_inverse_k = AI(A)
d_ik = d_ik_cubrix(2)

for n in range(0, 2**8):
    elements = [
            [[],
             []],

            [[],
             []]
    ]

    for k in range(1, 3):
        for i in range(1, 3):
            for j in range(1, 3):
                number = format(n, '#010b')[2:10]
                cell = float(number[4*(k-1) + 2*(i-1) + j-1])
                elements[k-1][i-1].append(cell)

    I = A_inverse_k
    J = Cubriz(elements)

    #(A.times(I, J) * 3).print()
    B = A.times(I, J) * 2

    if B.elements == d_ik.elements:
        print("Hurra")
"""

#(_A, D) = A.inverse_pair()
#_A.times(A, D).print()

D_ij = d_ij_cubrix(2)
D_jk = d_jk_cubrix(2)
D_ik = d_ik_cubrix(2)
d = D(2)

i = I(2)

#import pdb; pdb.set_trace()

#A.times(A**(-1), d).print()

for i in range(1, 3):
    for j in range(1, 3):
        for k in range(1, 3):
            txt = f"1 ="
            for l in range(1, 3):
                txt += f" + A_{i}{l}{k} I_{l}{j}{k} J_{i}{j}{l}"

            print(txt)

# Does it even make sense to ask whether there exist I and J such that IAJ = d_jk? Maybe not.
