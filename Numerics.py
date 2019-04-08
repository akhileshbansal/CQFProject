# -*- coding: utf-8 -*-

## module LUdecomp
import pandas as pd
import numpy as np


''' 
**** Code from Numerical Methods in Engineering with Python - by Jaan  Kiusalaas

    c,d,e = LUdecomp3(c,d,e).
    LU decomposition of tridiagonal matrix [c\d\e]. On output
    {c},{d} and {e} are the diagonals of the decomposed matrix
    x = LUsolve(c,d,e,b).
    Solves [c\d\e]{x} = {b}, where {c}, {d} and {e} are the
    vectors returned from LUdecomp3.
'''



def LUdecomp3(c,d,e):
    n = len(d)
    for k in range(1,n):
        lam = c[k-1]/d[k-1]
        d[k] = d[k] - lam*e[k-1]
        c[k-1] = lam
    return c,d,e

def LUsolve3(c,d,e,b):
    n = len(d)
    for k in range(1,n):
        b[k] = b[k] - c[k-1]*b[k-1]
    b[n-1] = b[n-1]/d[n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - e[k]*b[k+1])/d[k]
    return b


## module cubicSpline
''' k = curvatures(xData,yData).
Returns the curvatures of cubic spline at its knots.
y = evalSpline(xData,yData,k,x).
Evaluates cubic spline at x. The curvatures k can be
computed with the function ’curvatures’.
'''



def curvatures(xData,yData):
    n = len(xData) - 1
    c = np.zeros(n)
    d = np.ones(n+1)
    e = np.zeros(n)
    k = np.zeros(n+1)
    c[0:n-1] = xData[0:n-1] - xData[1:n]
    d[1:n] = 2.0*(xData[0:n-1] - xData[2:n+1])
    e[1:n] = xData[1:n] - xData[2:n+1]
    k[1:n] =6.0*(yData[0:n-1] - yData[1:n]) /(xData[0:n-1] - xData[1:n]) -6.0*(yData[1:n] - yData[2:n+1]) /(xData[1:n] - xData[2:n+1])   #ensure that K[0] = 0 and k[n] = 0 i.e extreme end points are having second diff as 0.0
    LUdecomp3(c,d,e)
    LUsolve3(c,d,e,k)
    return k


def evalSpline(xData,yData,k,x):
    
    def findSegment(xData,x):
        iLeft = 0
        iRight = len(xData)- 1
        while 1:
            if (iRight-iLeft) <= 1: return iLeft
            i =int((iLeft + iRight)/2)
            if x < xData[i]: iRight = i
            else: iLeft = i


    i = findSegment(xData,x)
    h = xData[i] - xData[i+1]
    y = ((x - xData[i+1])**3/h - (x - xData[i+1])*h)*k[i]/6.0 - ((x - xData[i])**3/h - (x - xData[i])*h)*k[i+1]/6.0 + (yData[i]*(x - xData[i+1]) - yData[i+1]*(x - xData[i]))/h
    return y

   



def jacobi(a, tol=1.0e-9):
    '''
        Jacobi method
     ****#http://w3mentor.com/learn/python/scientific-computation/
        #python-code-for-solving-eigenvalue-problem-by-jacobis-method/
    '''

    def maxElem(a):
        '''
            Find largest off-diag. element a[k,l]
        '''
        n = len(a)
        aMax = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(a[i, j]) >= aMax:
                    aMax = abs(a[i, j])
                    k = i
                    l = j
        return aMax, k, l

    def rotate(a, p, k, l):
        '''
            Rotate to make a[k,l] = 0
        '''
        n = len(a)
        aDiff = a[l, l] - a[k, k]
        if abs(a[k, l]) < abs(aDiff) * 1.0e-36:
            t = a[k, l] / aDiff
        else:
            phi = aDiff / (2.0 * a[k, l])
            t = 1.0 / (abs(phi) + np.sqrt(phi ** 2 + 1.0))
            if phi < 0.0:
                t = -t
        c = 1.0 / np.sqrt(t ** 2 + 1.0)
        s = t * c
        tau = s / (1.0 + c)
        temp = a[k, l]
        a[k, l] = 0.0
        a[k, k] = a[k, k] - t * temp
        a[l, l] = a[l, l] + t * temp
        for i in range(k):      # Case of i < k
            temp = a[i, k]
            a[i, k] = temp - s * (a[i, l] + tau * temp)
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(k + 1, l):  # Case of k < i < l
            temp = a[k, i]
            a[k, i] = temp - s * (a[i, l] + tau * a[k, i])
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(l + 1, n):  # Case of i > l
            temp = a[k, i]
            a[k, i] = temp - s * (a[l, i] + tau * temp)
            a[l, i] = a[l, i] + s * (temp - tau * a[l, i])
        for i in range(n):      # Update transformation matrix
            temp = p[i, k]
            p[i, k] = temp - s * (p[i, l] + tau * p[i, k])
            p[i, l] = p[i, l] + s * (temp - tau * p[i, l])

    n = len(a)
    maxRot = 5 * (n ** 2)    # Set limit on number of rotations
    p = np.identity(n) * 1.0    # Initialize transformation matrix
    for i in range(maxRot):  # Jacobi rotation loop
        aMax, k, l = maxElem(a)
        if aMax < tol:
            return np.diagonal(a), p
        rotate(a, p, k, l)
    print ('Jacobi method did not converge')
    
    
def calculate_pca(forwards, Period = 1, no_factors=3):
    fwddiff = forwards.diff(periods = Period)
    fwddiff = fwddiff.dropna()
    covmat = fwddiff.cov()
    covmat = covmat * (252/Period) / 10000
    eigenvecs, eigenmat = jacobi(covmat.values)
    eigvecs = pd.Series(eigenvecs, index=covmat.columns)
    sorted_eigvecs = eigvecs.sort_values(ascending=False)
    top3 = sorted_eigvecs[:no_factors].index
    eigenmat_df = pd.DataFrame(eigenmat, index=covmat.columns,
                            columns=covmat.columns)
    filtered_eigenmat = eigenmat_df.filter(top3)
    percent_variance = np.cumsum(sorted_eigvecs)/sum(sorted_eigvecs)
    return sorted_eigvecs, filtered_eigenmat, percent_variance
    
    
def trapInteg(FX, dt): #N points => N-1 segments to be added. This may be a bit less efficient than the other method but is just the one i use
    sumV = 0.0
    for i in range(1, len(FX)):
        sumV = sumV + (FX[i-1] + FX[i])
    sumV = sumV*dt*0.5
    return sumV 
