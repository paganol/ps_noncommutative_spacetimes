import numpy as np

#SS
def MSS(l,m):
    d = (2*l-1)*(2*l+3)
    first = np.sqrt(l*(l+1)/d)
    second = (l*l+l-3*m*m)/d
    return (first-second)*2/3*(-1)**m

def corrSS(cl,lmax):
    clout = np.zeros_like(cl)
    for l in range(2,lmax):
        for m in range(-l,l+1):
            clout[l] += MSS(l,m)*cl[l]
        clout[l] /= (2*l+1)
    return clout

def allcorrSS(cls,lmax):
    clsout = np.zeros_like(cls)
    for icl,cl in enumerate(cls):
        clsout[icl] = -3/16*corrSS(cl,lmax)
    return clsout
