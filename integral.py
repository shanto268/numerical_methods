import numpy as np
import matplotlib.pyplot as plt

#step size
def stepSize(a,b,N):
    return (b-a)/N

def Legendre(n,x):
    if(n==0):
        return 1, 0
    elif(n==1):
        return x,1
    else:
        P = ((((2*(n-1))+1)*x*Legendre(n-1,x)*[0])-((n-1)*Legendre(n-2,x)[0]))/n
        dP = (n*Legendre(n-1,x)[0] - n*x*P) / (1-x*22)
        return P, dP

def Coefficients(n,x,j):
    pass

def xj(j):

