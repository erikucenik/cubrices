import cubrices

print("""Welcome to the Cubrix playground!

We've provided a few tools for you to get a feel for the behavior
of cubrices. You can open this file inside a text editor and look
around for what might interest you.

If some Cubrix calls your attention, you can attach a .print(2) to
its end to see what values it holds.

If you still hold any doubts, don't forget you can always open
a Python interpreter and do:
    import cubrices
    help(cubrices)

Have fun!""")

# We'll be using cubic cubrices of size 3.
n = 2

# Cubrix generators
D_ij = cubrices.D_ij(n)
D_jk = cubrices.D_jk(n)
D_ik = cubrices.D_ik(n)
D    = cubrices.D(n)
I    = cubrices.I(n)
R    = cubrices.random(n, 0.0, 9.0)
F    = cubrices.Cubrix.from_function(n, lambda i, j, k: i + j - k)
A    = cubrices.Cubrix([
    [[1.0, 2.0],
     [3.0, 4.0]],

    [[5.0, 6.0],
     [7.0, 8.0]]
])

# Operator overloading
(D_ij + 3.5 + 1 + D_jk)
(F*3 + 3*F)
(I ** (-1))

# Printing
if False:
    # Rounded to 2
    str(I)
    print(I)

    # Rounded to 5
    I.to_str(5)
    I.print(5)

# Identity pairs. D is code for \Delta. Check the paper and help(cubrices.Cubrix).
# All of the following are equal to F.
(D_ij, I) = (cubrices.D_ij(n), cubrices.I(n))
(D_jk, I) = (cubrices.D_jk(n), cubrices.I(n))
(D_ik, I) = (cubrices.D_ik(n), cubrices.I(n))

F.times(D_ij, I)
F.times(I, D_jk)
D_ij.times(F, I)
I.times(F, D_ik)
D_jk.times(I, F)
I.times(D_ik, F)

F.times(D_ij, D_jk)
D_ij.times(F, D_ik)
D_jk.times(D_ik, F)

# Inverse pairs.
(A_inverse, D) = A.inverse_pair()

## i, j
A.times(A_inverse, D)
A_inverse.times(A, D)

## j, k
A.times(D, A_inverse)
A_inverse.times(D, A)

## i, k
D.times(A, A_inverse)
D.times(A_inverse, A)
