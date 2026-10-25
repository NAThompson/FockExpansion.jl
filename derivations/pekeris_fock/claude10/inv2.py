from monodromy import *
from mpmath import mpf, sqrt, nstr, matrix, eye, svd
ops = [[mp.mpf(x) for x in l.split()] for l in open('op_e2_A.txt')]
l0 = sqrt(mpf(8) / 9); sing = [l0, -l0, mpf(1), mpf(-1)]; base = mpf(1) / 5
T1 = monodromy_matrix(ops, base, l0, mpf(1) / 20, sing)
T2 = monodromy_matrix(ops, base, -l0, mpf(1) / 20, sing)
St = matrix(10, 5)
for i in range(5):
    for j in range(5):
        St[i, j] = (T1 - eye(5))[i, j]; St[i + 5, j] = (T2 - eye(5))[i, j]
print('singular values of [T1-1; T2-1]:', [nstr(x, 3) for x in svd(St, compute_uv=False)])
import pickle; pickle.dump({'T1': T1.tolist(), 'T2': T2.tolist()}, open('TA.pkl', 'wb'))
