import numpy as np

#lambda
def Mlambda(l,m):
    d = (2*l-1)*(2*l+3)
    first = np.sqrt(l*(l+1)/d)
    second = (l*l+l-3*m*m)/d
    third = (6*l-3*l*l+6*l**3+3*l**4+25*m*m-30*l*m*m-30*l*l*m*m+35*m**4)/d/(2*l-3)/(2*l+5)
    return ((first-second)*58/3+third*15/2)*(-1)**m

def corrlambda(cl,lmax):
    clout = np.zeros_like(cl)
    for l in range(2,lmax):
        for m in range(-l,l+1):
            clout[l] += Mlambda(l,m)*cl[l]
        clout[l] /= (2*l+1)
    return clout

def allcorrlambda(cls,lmax):
    clsout = np.zeros_like(cls)
    for icl,cl in enumerate(cls):
        clsout[icl] = corrlambda(cl,lmax)
    return clsout
