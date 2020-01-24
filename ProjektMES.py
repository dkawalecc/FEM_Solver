import numpy as np
import matplotlib.pyplot as plt
import math

inp = int(input("Podaj liczbę punktów: "))


def k(x, y):
    return 1


#    ____________
#   |4         3|
#   |    phi    |
#   |           |
#   |1_________2|
#
def phi(x, y, i):
    if i == 1:
        return (1 - x) * (1 - y)
    if i == 2:
        return x * (1 - y)
    if i == 3:
        return x * y
    if i == 4:
        return (1 - x) * y
    return 0


def derivativeX_of_phi(x, y, i):
    if phi(x, y, i) != 0:
        if i == 2 or i == 3:
            return 0.5
        if i == 1 or i == 4:
            return -0.5
    return 0


def derivativeY_of_phi(x, y, i):
    if phi(x, y, i) != 0:
        if i == 3 or i == 4:
            return 0.5
        if i == 1 or i == 2:
            return -0.5
    return 0


def e(x, y, i, f):
    # i - calculation node
    # a1 a2 are the width and height of element (surface)
    # b1 b2 are coords of the elements' starting point, bottom-left corner
    a1 = 1
    a2 = 1
    # E1
    if -1 <= x <= 0 and 0 <= y <= 1:
        b1 = -1
        b2 = 0
        if i == 0:
            return f((x - b1) / a1, (y - b2) / a2, 4)
        if i == 1:
            return f((x - b1) / a1, (y - b2) / a2, 1)
        return 0
    # E2
    if -1 <= x <= 0 and -1 <= y <= 0:
        b1 = -1
        b2 = -1
        if i == 1:
            return f((x - b1) / a1, (y - b2) / a2, 4)
        if i == 2:
            return f((x - b1) / a1, (y - b2) / a2, 1)
        if i == 3:
            return f((x - b1) / a1, (y - b2) / a2, 2)
        return 0
    # E3
    if 0 <= x <= 1 and -1 <= y <= 0:
        b1 = 0
        b2 = -1
        if i == 3:
            return f((x - b1) / a1, (y - b2) / a2, 1)
        if i == 4:
            return f((x - b1) / a1, (y - b2) / a2, 2)
        return 0
    return 0


def g(x, y):
    # r*sin(ϑ + π/2) = r*cos(ϑ) = x
    # y = rsin(ϑ)
    # g(x,y) = (x)^(2/3)
    return pow(x * x, 1 / 3)


def Bmatrix():
    B = [[0 for a in range(5)] for b in range(5)]
    # coords of the center points in certain areas: E1,E2,E3
    # variables a and b from e function => x = b1 + a1/2, y = b2 + a2/2
    coords = [(-0.5, 0.5), (-0.5, -0.5), (0.5, -0.5)]
    for i in range(5):
        for j in range(5):
            B[i][j] = sum([k(x, y) * e(x, y, i, derivativeX_of_phi) * e(x, y, j, derivativeX_of_phi)
                           + k(x, y) * e(x, y, i, derivativeY_of_phi) * e(x, y, j, derivativeY_of_phi)
                           for (x, y) in coords])
    return B


def Lmatrix():
    L = [0 for a in range(5)]
    # coords of the center points in certain areas: E1,E2,E3
    # Neumann clause(condition)
    coords = [(-0.5, 1), (-1, 0.5), (-1, -0.5), (-0.5, -1), (0.5, -1), (1, -0.5)]
    for i in range(5):
        L[i] = sum([k(x, y) * e(x, y, i, phi) * g(x, y) for (x, y) in coords])
    return L


def create_gauss_matrix(B, L, l):
    for curr_i in range(l - 1):
        max_row = max(range(curr_i, l - 1), key=lambda i: abs(B[i][curr_i]))
        (B[curr_i], B[max_row]) = (B[max_row], B[curr_i])
        (L[curr_i], L[max_row]) = (L[max_row], L[curr_i])
        for row in range(curr_i + 1, l):
            multiplier = B[row][curr_i] / B[curr_i][curr_i]
            B[row] = [0 for _ in range(curr_i + 1)] + [B[row][col] - multiplier * B[curr_i][col] for col
                                                       in range(curr_i + 1, l)]
            L[row] -= multiplier * L[curr_i]
    return (B, L)


def gauss(B, L):
    l = len(B)
    (B, L) = create_gauss_matrix(B, L, l)
    for curr_i in range(l - 1, 0, -1):
        for row_i in range(curr_i - 1, -1, -1):
            L[row_i] -= B[row_i][curr_i] / B[curr_i][curr_i] * L[curr_i]
            B[row_i][curr_i] = 0
    return [L[i] / B[i][i] for i in range(l)]


def solve():
    B = Bmatrix()
    L = Lmatrix()
    return gauss(B, L)


# u_h estimated
def u(w, x, y):
    return sum([w[i] * e(x, y, i, phi) for i in range(5)])


# inp - sample count
# 2 * inp so we can easily(aesthetically) exclude unnecessary surface from matrix
N = inp * 2 + 1
matrix = [[(0, 0) for _ in range(N)] for _ in range(N)]
inp_length = 1 / inp
w = solve()
print(w)

### do zmiany #####
for (x_start, y_start, i_start, j_start) in [(-1, 0, inp, 0), (-1, -1, 0, 0), (0, -1, 0, inp)]:
    y = y_start - inp_length / 2
    for i in range(i_start-1, i_start + inp):
        x = x_start + inp_length / 2
        for j in range(j_start, j_start + inp):
            matrix[i][j] = (x, y)
            x += inp_length
        y += inp_length
####################
T = [[0 for _ in range(N)] for _ in range(N)]
for (i, row) in enumerate(matrix):
    for (j, (x, y)) in enumerate(row):
        if not (x > 0 and y > 0):
            T[i][j] = u(w, x, y)

####################

X = np.linspace(-1, 1, N)
Y = np.linspace(-1, 1, N)

Xmesh, Ymesh = np.meshgrid(X, Y)

# T = np.random.randint(10, 30, (inp, inp))
# print(Ymesh)
# print(Xmesh)

plt.pcolormesh(X, Y, T, cmap='jet')
plt.colorbar()

plt.xlabel("X axis")
plt.ylabel("Y axis")

plt.show()
