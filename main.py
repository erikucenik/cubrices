import cubrices

# n is the size of the cubic cubrices.
n = 3

# Generates n-sized cubrices with random digits ranging from 0 to 9.
A = cubrices.random(n, 0.0, 9.0)
B = cubrices.random(n, 0.0, 9.0)

# Operator overloading
(A + B)
(3 * A)
(A ** (-1))  # Cubrix exponentiation (raises every element to a given power)

# Identity pairs. D is code for \Delta. Check the paper.
(D_ij, I) = (cubrices.D_ij(n), cubrices.I(n))
(D_jk, I) = (cubrices.D_jk(n), cubrices.I(n))

A.times(D_ij, I)
D_jk.times(I, A)

# Inverse pairs.
(A_inverse, D) = A.inverse_pair() # Returns (A^{-1}, \Delta)

A.times(A_inverse, D)
A.times(D, A_inverse)
C = D.times(A, A_inverse)

# You may print any of the previous cubrices with this method.
# The argument is the number of decimal places to round to.
C.print(3)
