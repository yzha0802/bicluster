import numpy as np
import scipy.linalg as la

def thred(z,delta):
    return np.sign(z)*(np.abs(z)>=delta)*(np.abs(z)-delta)

def ssvd(X,gamu = 2, gamv =2, merr = 10**(-4), niter = 100):
    n,d = X.shape
    #initial value of u and v
    U,s,VT = la.svd(X,full_matrices=False)
    u0 = U[:,0]
    v0 = VT.T[:,0]
    
    ud = 1
    vd = 1
    count = 0
    SST  = np.sum(X*X)
    
    while(ud>merr or vd>merr):
        count = count + 1
        
        # Update v
        z = X.T @ u0
        winv = np.abs(z)**gamv # weight inverse
        sigsq = np.abs(SST-np.sum(z*z))/(n*d-d)
        
        cand = z*winv  #candidate lambda
        delt = np.sort(np.append(np.abs(cand),0))
        delt_uniq = np.unique(delt)
        Bv = np.ones(len(delt_uniq)-1)*float("inf")
        
        ind = np.where(winv>10^(-8))
        cand1 = cand[ind]
        winv1 = winv[ind]
        for i in range(len(Bv)):
            vhat= thred(cand1,delta = delt_uniq[i])
            vhat = vhat/winv1
            vshrink = np.zeros(d)
            vshrink[ind] = vhat
            Bv[i] = np.sum((X - u0[:,None] @ vshrink[None,:])**2)/sigsq + np.sum(vhat!=0)*np.log(n*d)
        
        Iv = min(np.where(Bv== np.min(Bv))) #position of min BIC
        th = delt_uniq[Iv] #best lambda in this iteration
        vhat = thred(cand1,delta = th)
        vhat = vhat/winv1
        v1 = np.zeros(d)
        v1[ind] = vhat
        v1 = v1/(np.sqrt((np.sum(v1*v1)))) #v_new
        
        # Updating u
        z = X @ v1
        winu = np.abs(z)**gamu
        sigsq = np.abs(SST - np.sum(z*z))/(n*d-n)
        cand = z*winu
        delt = np.sort(np.append(np.abs(cand),0))
        delt_uniq = np.unique(delt)
        Bu = np.ones(len(delt_uniq)-1)*float("inf")
        ind = np.where(winu > 10^(-8))
        cand1 = cand[ind]
        winu1 = winu[ind]
        for i in range(len(Bu)):
            uhat = thred(cand1,delta = delt_uniq[i])
            uhat = uhat/winu1
            ushrink = np.zeros(n)
            ushrink[ind] = uhat
            Bu[i] = np.sum((X - ushrink[:,None] @ v1[None,:])**2)/sigsq + np.sum(uhat!=0)*np.log(n*d)
        Iu = min(np.where(Bu==np.min(Bu)))
        th = delt_uniq[Iu]
        uhat = thred(cand1,delta = th)
        uhat = uhat/winu1
        u1 = np.zeros(n)
        u1[ind] =  uhat
        u1 = u1/((np.sum(u1*u1))**0.5)
        
        
        ud = np.sqrt(np.sum((u0-u1)*(u0-u1)))
        vd = np.sqrt(np.sum((v0-v1)*(v0-v1)))
        
        if count > niter :
            print("Fail to converge! Increase the niter!")
            break
        
        u0 = u1
        v0 = v1
        
    s = u1[None, :] @ X @ v1[:, None] #ssvd layer is suv.T
    return u1, v1, s, count
