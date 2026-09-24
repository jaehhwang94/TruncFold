import math, numpy as np
from scipy.special import digamma, polygamma
from scipy.optimize import brentq

def estimate_scaled_f_prior(sample_vars,dfs):
    """Smyth-style scaled-F empirical-Bayes estimate from library variances."""
    s2=np.asarray(sample_vars,dtype=float); df=np.asarray(dfs,dtype=float)
    ok=np.isfinite(s2)&(s2>=0)&np.isfinite(df)&(df>0); s2=s2[ok]; df=df[ok]
    if not len(s2): raise ValueError('No variants with >=2 accepted poses; cannot estimate empirical-Bayes prior')
    s2=np.maximum(s2,np.finfo(float).tiny)
    z=np.log(s2)-digamma(df/2.0)+np.log(df/2.0)
    obs=float(np.var(z,ddof=1)) if len(z)>1 else 0.0
    evar=max(obs-float(np.mean(polygamma(1,df/2.0))),0.0)
    d0=1e8 if evar<=1e-12 else float(brentq(lambda d: float(polygamma(1,d/2.0)-evar),1e-4,1e8))
    log_s02=float(np.mean(z)+digamma(d0/2.0)-np.log(d0/2.0))
    return float(np.exp(log_s02)),d0

def truncfold_score(values,prior_var,prior_df):
    x=np.asarray(values,dtype=float); x=x[np.isfinite(x)]; n=len(x)
    if n==0: return None,None,None
    if n==1: return math.sqrt(prior_var),float('nan'),float(prior_var)
    s2=float(np.var(x,ddof=1)); d=n-1
    post=(prior_df*prior_var+d*s2)/(prior_df+d)
    return math.sqrt(post),s2,float(post)
