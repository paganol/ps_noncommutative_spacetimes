import numpy as np

#M1
def MM1_first(l,m):
    d = (2*l-1)*(2*l+3)
    first = np.sqrt(l*(l+1)/d)
    second = (l*l+l-3*m*m)/d
    return (first-second)*41/3*(-1)**m+11

def corrM1(cl,lmax):
    clout = np.zeros_like(cl)
    for l in range(2,lmax+1):
        for m in range(-l,l+1):
            clout[l] += MM1_first(l,m)*cl[l]
        clout[l] /= (2*l+1)
    return clout

def allcorrM1(cls,lmax):
    clsout = np.zeros_like(cls)
    for icl,cl in enumerate(cls):
        clsout[icl] = corrM1(cl,lmax)/32
    return clsout
