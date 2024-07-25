import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from spectraSS import allcorrSS
from spectraST import allcorrST
from spectraM1 import allcorrM1
from spectraM2 import allcorrM2
from spectralambda import allcorrlambda
from likelihood import exact_likelihood

doSS = False
doST = True
doM1 = False
doM2 = False
dolambda = False

fsky = 1
lmax = 2500

nas = 501
nrs = 1001

nrsover = 10001

incls = np.loadtxt('cls_tau_likelihood_0050.txt',usecols=(1,2,4)).T

clobs = incls


#SS
if doSS:
    As = np.linspace(0.995,1.005,nas)
    Rs = np.linspace(0,20,nrs)
    Rsover = np.linspace(0,20,nrsover)

    like = np.empty((nas,nrs))

    clcorr = allcorrSS(incls,lmax=lmax)

    for iA,A in enumerate(As):
        for iR,R in enumerate(Rs):
            clth = A*(incls+R*clcorr)
            like[iA,iR] = exact_likelihood(clobs,clth,fsky,lmax)

    plt.contour(Rs,As,like,[0,2.30,6.17])
    plt.show()

    prob2D = np.exp(-0.5*like)
    probpar = np.sum(prob2D,axis=0)*np.diff(As)[0]
    f = interp1d(Rs,probpar,kind='cubic')
    probover = f(Rsover)
    print(Rsover[np.where(probover.cumsum()/probover.sum()>0.9545)[0][0]])


#ST
if doST:
    As = np.linspace(0.95,1.05,nas)
    Rs = np.linspace(0,0.2,nrs)
    Rsover = np.linspace(0,0.2,nrsover)
    
    like = np.empty((nas,nrs))

    clcorr = allcorrST(incls,lmax=lmax)

    for iA,A in enumerate(As):
        for iR,R in enumerate(Rs):
            clth = A*(incls+R*clcorr)
            like[iA,iR] = exact_likelihood(clobs,clth,fsky,lmax)

    plt.contour(Rs,As,like,[0,2.30,6.17])
    plt.show()

    prob2D = np.exp(-0.5*like)
    probpar = np.sum(prob2D,axis=0)*np.diff(As)[0]
    f = interp1d(Rs,probpar,kind='cubic')
    probover = f(Rsover)
    print(Rsover[np.where(probover.cumsum()/probover.sum()>0.9545)[0][0]])


#M1
if doM1:
    As = np.linspace(0.995,1.005,nas)
    Rs = np.linspace(0,5,nrs)
    Rsover = np.linspace(0,5,nrsover)
    
    like = np.empty((nas,nrs))

    clcorr = allcorrM1(incls,lmax=lmax)

    for iA,A in enumerate(As):
        for iR,R in enumerate(Rs):
            clth = A*(incls+R*clcorr)
            like[iA,iR] = exact_likelihood(clobs,clth,fsky,lmax)

    plt.contour(Rs,As,like,[0,2.30,6.17])
    plt.show()

    prob2D = np.exp(-0.5*like)
    probpar = np.sum(prob2D,axis=0)*np.diff(As)[0]
    f = interp1d(Rs,probpar,kind='cubic')
    probover = f(Rsover)
    print(Rsover[np.where(probover.cumsum()/probover.sum()>0.9545)[0][0]])


#M2
if doM2:
    As = np.linspace(0.995,1.005,nas)
    Rs = np.linspace(0,10,nrs)
    Rsover = np.linspace(0,10,nrsover)

    like = np.empty((nas,nrs))

    clcorr = allcorrM2(incls,lmax=lmax)

    for iA,A in enumerate(As):
        for iR,R in enumerate(Rs):
            clth = A*(incls+R*clcorr)
            like[iA,iR] = exact_likelihood(clobs,clth,fsky,lmax)

    plt.contour(Rs,As,like,[0,2.30,6.17])
    plt.show()

    prob2D = np.exp(-0.5*like)
    probpar = np.sum(prob2D,axis=0)*np.diff(As)[0]
    f = interp1d(Rs,probpar,kind='cubic')
    probover = f(Rsover)
    print(Rsover[np.where(probover.cumsum()/probover.sum()>0.9545)[0][0]])

#lambda
if dolambda:
    As = np.linspace(0.995,1.005,nas)
    Rs = np.linspace(0,0.1,nrs)
    Rsover = np.linspace(0,0.1,nrsover)
    
    like = np.empty((nas,nrs))

    clcorr = allcorrlambda(incls,lmax=lmax)

    for iA,A in enumerate(As):
        for iR,R in enumerate(Rs):
            clth = A*(incls+R*clcorr)
            like[iA,iR] = exact_likelihood(clobs,clth,fsky,lmax)

    plt.contour(Rs,As,like,[0,2.30,6.17])
    plt.show()

    prob2D = np.exp(-0.5*like)
    probpar = np.sum(prob2D,axis=0)*np.diff(As)[0]
    f = interp1d(Rs,probpar,kind='cubic')
    probover = f(Rsover)
    print(Rsover[np.where(probover.cumsum()/probover.sum()>0.9545)[0][0]])

