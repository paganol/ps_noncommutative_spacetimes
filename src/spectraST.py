import numpy as np

#ST
def MST(l,m):
    d = (2*l-1)*(2*l+3)
    first = 0.5*np.sqrt(l*(l+1)/d)
    second = (l*l+l-3*m*m)/d
    return (first+second)*58/3*(-1)**m + 11

def corrST(cl,lmax):
    clout = np.zeros_like(cl)
    for l in range(2,lmax+1):
        for m in range(-l,l+1):
            clout[l] += MST(l,m)*cl[l]
        clout[l] /= (2*l+1)
    return clout

def allcorrST(cls,lmax):
    clsout = np.zeros_like(cls)
    for icl,cl in enumerate(cls):
        clsout[icl] = corrST(cl,lmax)/32
    return clsout
