#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 4 14:59:15 2020

@author: cwuyichen
clustering analysis main program
following CCToolbox(http://www.datalab.uci.edu/resources/CCT/)
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statsmodels.api as sm

def wls(y,x,w):
    # y: 1D;  
    # w: 1D;
    # x: array;
    modelResult = sm.WLS(y,x,w).fit()
    mu = modelResult.params
    res = modelResult.resid
    sigma = np.dot(w*res,res)/w.sum()
    return mu,sigma

def get_pdf(track,reg_track,sigma): 
    d = track.shape[1]
    isigma = np.linalg.inv(sigma)
    dsigma = np.linalg.det(sigma)
    
    track = track-reg_track
    ans = np.exp(-0.5*np.sum(np.dot(track,isigma)*track,axis=1)) 
    ans = ans*np.power(2*np.pi,-d/2.)*np.power(dsigma,-1/2.)
    return ans   # [tcnum,7]

def get_w(p,dots,alpha):
    ans = np.zeros((len(dots),len(alpha)))
    st = 0
    for i in range(len(dots)):
        ans[i,:] = np.prod(p[st:st+dots[i],:],axis=0)*alpha
        st = st+dots[i]
    return ans  # [tcnum,7]


def findNan(arr,a):
    sumarr = np.sum(arr,axis=1)
    lis = np.argwhere(sumarr==0.0)
    if len(lis) != 0:
        for i in lis:
            arr[i] = np.finfo(float).tiny*1e100*a   
            
def weightNorm(w):
    dummy = np.ones(w.shape)
    dummy = np.multiply(dummy.T,w.sum(axis=1))
    return np.divide(w,dummy.T)   

def stop_check(l,itern,limit=50,stopval=np.power(0.1,8)):
    if itern>=limit:
        #print(f'{itern}th iterration. Limited iterration.')
        return False
    if itern>1:
        if np.isnan(l[itern-1]):
            #print(f'{itern}th iterration. Likelihood is nan.')
            return False
        if l[itern-1]<l[itern-2]:
            #print(f'{itern}th iterration. Likelihood is descresing.')
            return False
        else:
            abs_change = l[itern-1]-l[0]
            if abs_change == 0:
                #print(f'{itern}th iterration. Likelihood is equal to the first iteration.')
                return False
            else:
                delta = l[itern-1]-l[itern-2]
                if delta/abs_change < stopval:
                    #print(f'{itern}th iterration. Likelihood is converged.')
                    return False
                else:
                    #print(f'{itern}th iterration.')
                    return True
    else:
        #print('First iterration'.)
        return True
        
#%% data preparing------------------------------------------------------------
p = 2   # order of polynomial regression
c = 6   # cluster number
iternum = 50
emnum = 10
filename = './trackData.csv'

with open(filename,'r') as f:
    df = pd.read_csv(f)
      
    dots_oftc = np.array(df.groupby('number').size())
    tcnum = len(dots_oftc)  
    dotnum = len(df)  

# lat and lon  [dotnum,2]
Y = np.array(df[['lon','lat']])  
  
# regression matrices  [dotnum,3]
hasX = False
for i in dots_oftc:
    t = np.vander(np.arange(i), N=p+1,increasing=True)
    if hasX:
        X = np.concatenate((X,t),axis=0)
    else:
        X = t
        hasX = True
print(f'There are {len(dots_oftc)} TCs.')


#%%  regression ----------------------------------------------------   
final_c = []
final_likelihood = -np.inf
final_likelihoodp = 0
for em in range(emnum):
    print('emstart {}'.format(em))
    # weights  [dotnum,7]             
    # np.random.seed(10)
    weight = -0.5*np.log(np.random.rand(tcnum,c))    # just follow CCToolbox 
    weight = weightNorm(weight)

    try:
        go = True
        iternow = 0
        likelihood = np.zeros(iternum)
        while go: 
            iternow = iternow+1
            # M-step: re-calculate the regression matrix (B and S)------------------
            W = np.zeros((dotnum,c))
            for k in range(c):
                st = 0
                for i in range(tcnum):
                    W[st:st+dots_oftc[i],k] = weight[i,k]
                    st = st+dots_oftc[i]
            A = np.sum(weight,axis=0)/tcnum
            
            B = np.zeros((p+1,Y.shape[1],c))
            S = np.zeros((Y.shape[1],Y.shape[1],c))
            for k in range(c):
                B[:,0,k],S[0,0,k] = wls(Y[:,0],X,W[:,k])   # 0 for longitude
                B[:,1,k],S[1,1,k] = wls(Y[:,1],X,W[:,k])   # 1 for latitude
            
            # E-step: compute the membership (weight) ---------------------------
            pdf = np.zeros(W.shape)
            for k in range(c):
                pdf[:,k] = get_pdf(Y,np.dot(X,B[:,:,k]),S[:,:,k])
            scale = pdf.mean()
            pdf = pdf/scale           # normalize the pdf with a view to its mean
            weight = get_w(pdf,dots_oftc,A)
            
            # check the stop conditions ------------------------------------------
            findNan(weight,A)     # find the zero in weight and make it the possible minimun value in order to avoid error
            weight_sum = weight.sum(axis=1)
            likelihood[iternow-1] = np.log(weight_sum).sum()+pdf.shape[0]*np.log(scale)     # add log(scale) to diminish the effect of normalization
            weight = weightNorm(weight)
            
            go = stop_check(likelihood,iternow,iternum)
    except:
        print('error')
        continue
    result_c = np.argmax(weight,axis=1)+1     
    result_likelihood = likelihood[iternow-1]
    result_likelihoodp = result_likelihood/Y.shape[0]/Y.shape[1]
    print(f'emstart {em}, likelihood={result_likelihood:.2f}')
    
    if result_likelihood > final_likelihood:
        final_likelihood = result_likelihood
        final_c = result_c
        final_likelihoodp = result_likelihoodp
        print(f'{em} emstart is chosen.')
        
print(filename)
print(f'cluster number: {c}, likelihood:{final_likelihoodp:.4f}')
print('remember to save the result!')

#%% output ----------------------------------------------
outf = './clusterResult.csv'.format(c)
np.savetxt(outf, final_c,fmt="%2d", delimiter=",")

