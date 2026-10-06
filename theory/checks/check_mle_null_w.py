"""Check V1: does the Gaussian MLE of w fail when some w = 0?

Reproduces the design of v6 slide 136 (Scheme 2.1): 20 SNPs in two blocks of 10,
q = 200 traits, N = 10,000, h2 = 0.8, 10% of SNPs with w = 0 and the rest scaled
F(100, 100). Effects are Gaussian given w, so z is exactly multivariate normal.
The slides reported that the MLE fails here (R2 0.068, null SNPs at 0.5-1.6).
With an exact likelihood, analytic gradient and L-BFGS-B bounds w >= 0 the MLE
puts every null SNP at 0 and beats OLS, so the slide's failure was numerical,
not a consequence of non-normality. LD is AR(1) with r = 0.6 within each block.

Run: python3 check_mle_null_w.py
"""
import numpy as np
from scipy.optimize import minimize
rng=np.random.default_rng(1)
p,q,N,h2=20,200,10000,0.8
def ar1(k,r): i=np.arange(k); return r**np.abs(i[:,None]-i[None,:])
R=np.zeros((p,p)); R[:10,:10]=ar1(10,0.6); R[10:,10:]=ar1(10,0.6)
s=N*h2/p
def sim():
    w=rng.f(100,100,p); w=w/ (100/98); z0=rng.random(p)<0.1; w[z0]=0; w=w/w.mean()
    B=rng.normal(size=(p,q))*np.sqrt(w*h2/p)[:,None]
    L=np.linalg.cholesky(R)
    Z=np.sqrt(N)*(R@B)+L@rng.normal(size=(p,q))
    return w,Z
def nll(w,Z):
    S=R+s*R@np.diag(w)@R
    sign,ld=np.linalg.slogdet(S); Si=np.linalg.inv(S)
    return 0.5*(q*ld+np.einsum('ia,ij,ja->',Z,Si,Z))
def grad(w,Z):
    S=R+s*R@np.diag(w)@R; Si=np.linalg.inv(S)
    A=Si@Z; M=q*Si-A@A.T   # d nll/dS = 0.5*(q Si - Si ZZ' Si)
    RMR=R@M@R
    return 0.5*s*np.diag(RMR)
D=R*R
out=[]
for rep in range(50):
    w,Z=sim()
    res=minimize(nll,np.ones(p),args=(Z,),jac=grad,method='L-BFGS-B',bounds=[(0,None)]*p)
    ols=np.linalg.solve(D,(Z**2-1).sum(1)*s)/(q*s*s)
    nul=w==0
    out.append((res.success,res.x[nul].mean() if nul.any() else np.nan,np.corrcoef(w,res.x)[0,1]**2,np.corrcoef(w,ols)[0,1]**2,(res.x[nul]<0.05).mean() if nul.any() else np.nan))
o=np.array(out,dtype=float)
print('reps',len(o),'converged',int(o[:,0].sum()))
print('mean MLE w-hat at true w=0: %.3f; share of nulls with MLE<0.05: %.2f'%(np.nanmean(o[:,1]),np.nanmean(o[:,4])))
print('mean R2 (w, w-hat): MLE %.3f, OLS %.3f'%(o[:,2].mean(),o[:,3].mean()))
