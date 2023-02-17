import numpy as np
import time
from scipy.linalg import lu, svd


np.random.seed(10)

# ランダムな要素をもつ行列を生成
n = 10001
A = np.random.randint(1, 10, (n, n))

# PA=LU分解
c = time.time()
P, L, U = lu(A)
print(time.time() - c)

c = time.time()
U, s, V = np.linalg.svd(A, full_matrices=True)
print(time.time() - c)

# print("A:\n{}\n".format(A))
# print("P:\n{}\n".format(P))
# print("L:\n{}\n".format(L))
# print("U:\n{}".format(U))

from OCC.Core.math import math_Matrix, math_Vector
from OCC.Core.math import math_GaussLeastSquare, math_SVD

vec_b = math_Vector(1, n)
for i in range(n):
    vec_b.SetValue(i + 1, 1.0)
vec_x = math_Vector(1, n)
for i in range(n):
    vec_x.SetValue(i + 1, 0.0)
mat = math_Matrix(1, n, 1, n)
for (i, j), v in np.ndenumerate(A):
    mat.SetValue(i + 1, j + 1, float(v))

c = time.time()
mat_LU = math_GaussLeastSquare(mat)
mat_LU.Solve(vec_b, vec_x)
print(time.time() - c)

c = time.time()
mat_SVD = math_SVD(mat)
mat_SVD.Solve(vec_b, vec_x)
print(time.time() - c)

# print(vec_b)
# vec_txt = ""
# for i in range(vec_b.Lower(), vec_b.Upper() + 1):
#    vec_txt += f"{vec_b.Value(i)}\t"
# print(vec_txt)
# print(vec_x)
# vec_txt = ""
# for i in range(vec_x.Lower(), vec_x.Upper() + 1):
#    vec_txt += f"{vec_x.Value(i)}\t"
# print(vec_txt)
#
