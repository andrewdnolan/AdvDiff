#!/usr/local/bin/python3

"""
Iterative Solver Functions
written as part of MATH 709

author: andrewnolan

To Do:
  [ ]  sparse matrix implementation
  [ ]  add conjugate gradient
  [ ]  add GMRES
"""

import numpy as np
from scipy import linalg as LA

def jacobi(A,b):
    L = np.tril(-A,-1)
    D = np.diag(np.diag(A))
    U = np.triu(-A,1)

    T = LA.inv(D) @ (L+U)
    C = LA.inv(D) @ b

    itteration_limit = 10000

    x = np.zeros_like(b)
    x_new = np.zeros_like(b)

    for __ in range(itteration_limit):
        x_new = T @ x + C
        if np.allclose(x,x_new,rtol=1e-10):
            break
        x = x_new

    return x_new

def gauss_siedel(A,b):
    L = np.tril(-A,-1)
    D = np.diag(np.diag(A))
    U = np.triu(-A,1)

    T = LA.inv(D-L) @ U
    C = np.dot(LA.inv(D-L),b)

    x = np.zeros_like(b)
    x_new = np.zeros_like(b)

    itteration_limit = 10000

    for __ in range(itteration_limit):
        x_new = T @ x + C
        if np.allclose(x,x_new,rtol=1e-10):
                break
        x = x_new
    return x_new

def SOR(A,b,omega=1.75):
    L = np.tril(-A,-1)
    D = np.diag(np.diag(A))
    U = np.triu(-A,1)


    T = LA.inv(D-omega*L)@(omega*U+(1-omega)*D)
    C = omega*LA.inv(D-(omega*L)) @ b

    x = np.zeros_like(b)
    x_new = np.zeros_like(b)

    itteration_limit = 10000

    for __ in range(itteration_limit):
        x_new = T@x + C
        if np.allclose(x,x_new,rtol=1e-10):
            break
        x = x_new

    return x_new

def TDMA(A,d):
    """Tridiagonal matrix algorithm (Thomas Algorithm)
    ref:https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    """
    a  = np.diag(A, -1).copy()
    b  = np.diag(A).copy()
    c  = np.diag(A, 1).copy()
    d  = d.copy()
    ne = len(d) # number of equations

    for k in range(1,ne):
        mm   = a[k-1] / b[k-1]
        b[k] = b[k] - mm * c[k-1]
        d[k] = d[k] - mm * d[k-1]

    x = d
    x[-1] = d[-1]/b[-1]

    for k in range(ne-2,-1,-1):
        x[k] = (d[k] - c[k]*x[k+1])/b[k]

    return x
