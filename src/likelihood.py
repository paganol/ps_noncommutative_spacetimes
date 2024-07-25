import numpy as np

def exact_likelihood(clobs,clth,fsky,lmax):

    assert clobs.shape == clth.shape
    like = 0
    if clobs.ndim == 1:
        assert len(clobs) > lmax
        l = np.arange(lmax+1)
        like = fsky*np.sum((2*l[2:lmax+1]+1)*(clobs[2:lmax+1]/clth[2:lmax+1]+np.log(clth[2:lmax+1])-np.log(clobs[2:lmax+1])-1))
    else:
        assert len(clobs[0]) > lmax
        l = np.arange(lmax+1)
        detth = clth[0,2:lmax+1]*clth[1,2:lmax+1]-clth[2,2:lmax+1]*clth[2,2:lmax+1]
        detobs= clobs[0,2:lmax+1]*clobs[1,2:lmax+1]-clobs[2,2:lmax+1]*clobs[2,2:lmax+1]
        tr = (clobs[0,2:lmax+1]*clth[1,2:lmax+1]+clobs[1,2:lmax+1]*clth[0,2:lmax+1]-2*clobs[2,2:lmax+1]*clth[2,2:lmax+1])/detth
        like = fsky*np.sum((2*l[2:lmax+1]+1)*(tr+np.log(detth)-np.log(detobs)-2))
    return like
