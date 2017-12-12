import numpy as np
import matplotlib.pyplot as p
from os import system,path
from joblib import Parallel,delayed
import time
#import modred as mr
from scipy.interpolate import rbf
from scipy.stats import linregress
from scipy.optimize import fmin,differential_evolution,fmin_cg,fmin_powell
import sys
import time
# sys.path.append('../../Utilities/')
# import surrogate as sg
# import LHD as lhd
# import PSO as pso
# reload(sg)
# reload(lhd)
# reload(pso)
#sys.path.append('../../py_rx/')
if path.exists('pythonrxNOdx.so'):
  print "****Remove module"
  system('rm pythonrxNOdx.so')
system('f2py --fcompiler=gfortran -c ./umat_1d_nograins.f -m pythonrxNOdx')
if path.exists('pythonrxNOdx.so'):
  print " * * * * * * * * COMPILE SUCCESS"
import pythonrxNOdx as rx
reload(rx)

p.ion()
p.close('all')

onlyharden,dowhat =True,'rx2dx_hard'
onlyharden,dowhat =False,'rx2dx_all'

cutfigs = False#True
minmax = False#True
absdiff = True

usetests = True
fvmaxstresspt = False
nrpypts = 100.

fnum=1

doplots= True#False
doplotsfinal =True
do_opt = False
printopt = True
useNelder = True#False
pydir = './python/'
#pydir = './abq_nl/'
figdir = pydir+'figures/'+dowhat+'/'
if path.exists(pydir+'figures/')==False:
  system('mkdir '+pydir+'figures/')
if path.exists(figdir)==False:
  system('mkdir '+figdir)
 
uselogstrain = False #linear case =False
usetannermod = False

## ALL THE CURVES:
crate = [0.0004]*10 + [6000.,1.,0.1,0.01] + [5200.,1.,0.1,0.01] + [1.,0.01]
ctemp = [25.,134., 202., 236., 269., 286., 303., 337., 405., 541.] + [25.]*4 + [269.]*4 + [541.]*2
ccolor = ['k','r','b','g','c','m','y']*10
ccolor = ['k']*30
csymb = ['-']*27+['--']*7+['-.']*7+[':']*7

# ## WITHOUT HIGH RATE:
crate = [0.0004]*10 + [1.,0.1,0.01] + [1.,0.1,0.01] + [1.,0.01]
ctemp = [25.,134., 202., 236., 269., 286., 303., 337., 405., 541.] + [25.]*3 + [269.]*3 + [541.]*2
ccolor = ['k','r','b','g','c','m','y']*10
ccolor = ['k']*30
csymb = ['-']*27+['--']*7+['-.']*7+[':']*7

## without high temp
#crate = [0.0004]*9# + [1.,0.1,0.01] + [1.,0.1,0.01]# + [1.,0.01]
#ctemp = [25.,134., 202., 236., 269., 286., 303., 337., 405.]# + [25.]*3 + [269.]*3 #+ [541.]*2



### ONLY 3 temps
# crate = [0.0004]*3 + [6000.,1.,0.1,0.01] + [5200.,1.,0.1,0.01] + [1.,0.01]
# ctemp = [25.,269.,541.] + [25.]*4 + [269.]*4 + [541.]*2
# ccolor = ['k','r','b']+['k']*4+['r']*4+['b']*2
# csymb = ['']*27+['--']*7+['-.']*7+[':']*7

if onlyharden:
  crate = [1.,0.1,0.01,0.0004,0.0004,5200.,1.,0.1,0.01]
  ctemp = [25.,25.,25.,25.,134.,269.,269.,269.,269.]
  ccolor = ['k']*4 + ['b'] + ['r']*4
  
  crate = [1.,0.1,0.01,0.0004,0.0004,0.0004,1.,0.1,0.01]
  ctemp = [25.,25.,25.,25.,134.,202.,269.,269.,269.]
  ccolor = ['k']*4 + ['b'] + ['g']+ ['r']*3
  
  
#only a few:
#crate = [0.0004]*8 + [1.,0.1,0.01] + [1.,0.1,0.01] + [1.]
#ctemp = [25.,134., 202., 236., 269., 286., 303., 337.] + [25.]*3 + [269.]*3 + [541.]

nrdata = crate.__len__()
cnames = ["curve_%s_%s"%(ctemp[i]+273,crate[i]) for i in range(nrdata)]
nra,nte = 'rate','temp'
testnames = ['542_0.0004to0.1','298_6000to0.0004','298_0.0004to542_1min','298_0.0004to542_20min','542_0.0004to298','298_6000to542_0.0004_45sec','542_5200to298_0.0004' ]
testnames = ['542_0.0004to0.1','298_0.0004to542_1min','298_0.0004to542_20min','542_0.0004to298']
#testnames = ['542_0.0004to0.1','542_0.0004to298']

nrtests = testnames.__len__()

data = np.ma.load('data_Tanner') # order = np.c_[times,strains,stresses,temps,rates]
testdata = np.ma.load('testdata_Tanner') # order = np.c_[times,strains,stresses,temps,rates,stepnumber]

Tk0 = np.array([25,269,405,541,691]*2)+273.
muD0 = np.array([42000,38110,35940,33460,31370,42130,38080,35820,33560,31070])

Tk = np.array(ctemp)+273.
rates = np.array(crate)


#
#
#
#
#
#
#
#
#  PAPER VALUES
#
#
#
#
V = {'a0': 2.1037,
    'a02': 0.904,
    'a03': 6370.675,
     'c0': 0.,
     'c1': 278.87,
    'c20': 11.773,
    'c30': 83.07,
     'c4': 6378.74,
    'c50': 9346.62,
    'c5c': 393.44,
   'c5g0': 47.247,
   'c5t0': 17634,
    'dx0': 1,
   'dx02': 1,
    'dxh': 1,
  'dxphi': 1,
    'dxu': 1,
    'ed0': 4.7e3,
   'emu0': 43.8e3,
    'enu': 1./3,
    'et0': 252.,
      'p': 1.,
      'q': 2.,
     'r3': 0.8079,
     'r4': 0.769,
    'r5a': 0.1424,
    'r5b': 1.7677,
    'r5g': 3.87,
  'rate0': 1e6,
 'rate02': 4.0112e12,
   'sig0': 17.295,
   'siga': 12.519}
   

keynames = ["emu0","ed0","et0","enu","siga","sig0","dx0","dxu",'dxh',"dxphi",
        "c0","c1","c20","a02","rate02","c30","r3","a03","c4","r4",
        "c50","c5t0","c5g0","r5g","r5a","r5b","c5c",
        "a0","rate0","q","p","dx02"]


nrstatv=500
dx0 = 300.

dtimemax0 = 10
dtimemaxtest = 10
dtimemax = dtimemax0

def getPYstress(strains0,times0,temps0,**V):

  time_ref0 = time.time()
  props = np.abs(np.array([V[kn] for kn in keynames]))
  totalresp = np.array([])
  
  PYstrain,PYstress,PYdx,PYeqpl = [],[],[],[]
  # initialise state variables:
  statev = np.zeros(nrstatv)
#  statev[4] = dx0
#  statev[5] = dx0
  statev[4] = 1.
#  statev[9] = dx0
# statev[10] = dx0
  statev[6]=0.01
  # initial total strain
  totalstrain = 0.
  stress,stressarr = 0.,[0.]
  for i in range(strains0.size-1):
    dstrain = strains0[i+1]-strains0[i]
    dtime = times0[i+1]-times0[i]
    temp = temps0[i+1] #(temps0[i+1]+temps0[i])/2.
    dtemp = (temps0[i+1]-temps0[i])/2.
    nsteps = 1
    if dtime>dtimemax:
      nsteps = np.ceil(dtime/dtimemax)
      #print nsteps
      dtime /= nsteps
      dtemp /= nsteps
      dstrain /= nsteps
      nsteps = int(nsteps)
    for nstp in range(nsteps):
      stress,statev = rx.umat1d(stress,statev.copy(),dstrain,dtime,temp,props)
    stressarr += [stress]
  return np.array(stressarr)
     


Xnames = ['a0','a02','c1','c20','sig0','c4']


Xnames = ['sig0','a02','c1','c20','c4']#,'r4']

Xnames += ['ed0','emu0']
#Xnames += ["c50",'c5t0',"c5g0","c5c",'r5a','r5b','r5g']

Xnames += ['c30','a03']

Xnames += ['p','q','r4']
Xnames += ['a0','rate02']

#Xnames += ['siga']

Xnames += ["c50",'c5t0',"c5g0"]
Xnames += ["c5c",'r5a','r5b','r5g']

X0 = np.array([V[Xnames[xnr]] for xnr in range(Xnames.__len__())])
SF = 10**(np.log10(X0+10))
  

Xbounds = []
for xnr in range(X0.size):
  Xbounds += [(X0[xnr]/2,X0[xnr]*2)]
  
X0 /= SF

consf = []

cABS = lambda x: np.prod(x[np.where(x<>0.)[0]]/np.abs(x[np.where(x<>0.)[0]]) > 0) # all values have to be positive
consf = [cABS]

mu = lambda T: V['emu0'] - V['ed0']/(np.exp(V['et0']/T)-1.)
#
dopval = False
doqval = False
if 'p' in Xnames:
 dopval = True
 pnr = np.where([Xnames[i]=='p' for i in range(Xnames.__len__())])[0][0]
 cp = lambda x: np.abs(x[pnr]*SF[pnr] - 0.5) <= 0.5 # 0<p<1
 Xbounds[pnr] = (0.,1.)
 consf += [cp]
 
if 'q' in Xnames:
 doqval = True
 qnr = np.where([Xnames[i]=='q' for i in range(Xnames.__len__())])[0][0]
 cq = lambda x: np.abs(x[qnr]*SF[qnr]-1.5) <= 0.5 # 1<q<2
 Xbounds[qnr] = (1.,2.)
 consf += [cq]
 
if 'a0' in Xnames:
  a0nr = np.where([Xnames[i]=='a0' for i in range(Xnames.__len__())])[0][0]
  if dopval*doqval:
    sfe = lambda x: (1.-((200./(x[a0nr]*SF[a0nr]*mu(200.)))*np.log(V['rate0']/1.e-8))**(1./(SF[qnr]*x[qnr])))**(1./(x[pnr]*SF[pnr]))<1.
  elif dopval:
    sfe = lambda x: (1.-((200./(x[a0nr]*SF[a0nr]*mu(200.)))*np.log(V['rate0']/1.e-8))**(1./V['q']))**(1./(SF[pnr]*x[pnr]))<1.
  elif doqval:
    sfe = lambda x: (1.-((200./(SF[a0nr]*x[a0nr]*mu(200.)))*np.log(V['rate0']/1.e-8))**(1./(SF[qnr]*x[qnr])))**(1./V['p'])<1.
  else:
   sfe = lambda x: (1.-((200./(SF[a0nr]*x[a0nr]*mu(200.)))*np.log(V['rate0']/1.e-8))**(1./V['q']))**(1./V['p'])<1.
  consf += [sfe]

if 'a02' in Xnames:
  a02nr =  np.where([Xnames[i]=='a02' for i in range(Xnames.__len__())])[0][0]
  c2 = lambda x: np.exp((1700./(SF[a02nr]*x[a02nr]*mu(1700.)))*np.log(V['rate02'])/1.e-8)>1.
  consf += [c2]

if "r4" in Xnames:
  r4nr =  np.where([Xnames[i]=='r4' for i in range(Xnames.__len__())])[0][0]
  cr4 = lambda x: np.abs(SF[r4nr]*x[r4nr]-0.75)<=0.25
  Xbounds[r4nr] = (0.5,1.)
  consf += [cr4]
  

if "r5a" in Xnames:
  r5anr =  np.where([Xnames[i]=='r5a' for i in range(Xnames.__len__())])[0][0]
  cr5a = lambda x: SF[r5anr]*x[r5anr]<=1.
  Xbounds[r5anr] = (0.,1.)
  consf += [cr5a]
  
for xnr in range(X0.size):
  Xbounds[xnr] /= SF[xnr]#[(X0[xnr]/2,X0[xnr]*2)]
  
# if "r5b" in Xnames:
#   r5bnr =  np.where([Xnames[i]=='r5b' for i in range(Xnames.__len__())])[0][0]
#   cr5b = lambda x: x[r5bnr]>=1.
#   consf += [cr5b]

     
 
# if doplots:
#   p.figure(1)
#   for cnr in range(nrdata):# [1]:
#     cnc = cnames[cnr]
#     curve = data[cnc]
#     times0 = curve[:,0]
#     strains0 = curve[:,1]
#     stresses0 = curve[:,2]
#     temps0 = curve[:,3]
#     pystress = getPYstress(strains0,times0,temps0,**V)
#     p.plot(strains0,stresses0,'k.',ms=2)#ccolor[cnr]+'.')
#     p.plot(strains0,pystress,'k-')#ccolor[cnr])
#  
#   p.xlabel(r'$\varepsilon$',fontsize=16.)    
#   p.ylabel(r'$\sigma$',fontsize=16.)
# #p.title('OFHC Cu Data')
# p.title('Optimised Parameter Fit')
# p.legend(['Data','Model'],loc='lower right')
 
 


def objF(X):
  Fv=0.
  Xeff = np.abs(X)*SF
  for xnr in range(X.size):
    V[Xnames[xnr]]=Xeff[xnr]#np.abs(X[xnr])*SF[xnr]
  
  if np.prod([consf[i](X) for i in range(consf.__len__())]) == 0:
    Fv = 100000.
  else:
    dtimemax = dtimemax0
    for cnr in range(nrdata):# [1]:
      cnc = cnames[cnr]
      curve = data[cnc]
      times0 = curve[:,0]
      strains0 = curve[:,1]
      stresses0 = curve[:,2]
      temps0 = curve[:,3]
      if cutfigs:
        mxpt = np.where(stresses0==np.max(stresses0))[0][0]
        strains0 = strains0[:mxpt]
        stresses0 = stresses0[:mxpt]
      pystress = getPYstress(strains0,times0,temps0,**V)
      if absdiff:
	diffs = (pystress[1:]-stresses0[1:])#/stresses0[1:]
      else:
	diffs = (pystress[1:]-stresses0[1:])/stresses0[1:]
      #print cnr
      #print pystress[1:]
      #print stresses0[1:]
      if fvmaxstresspt:
        str0 = strains0[np.where(stresses0==np.max(stresses0))[0]]
        str1 = strains0[np.where(pystress==np.max(pystress))[0]]
	Fv += 100*np.abs(str0-str1)
      elif minmax:
        Fv = np.max([Fv,1000*np.max(np.abs(diffs))])
      else:
        Fv += np.sum(diffs*diffs)/diffs.size
      
    if usetests:
      dtimemax = dtimemaxtest
      for cnr in range(nrtests):# [1]:
        cnc = testnames[cnr]
        curve = testdata[cnc]
        times0 = curve[:,0]
        strains0 = curve[:,1]
        stresses0 = curve[:,2]
        temps0 = curve[:,3]
        pystress = getPYstress(strains0,times0,temps0,**V)
	if absdiff:
	  diffs = (pystress[1:]-stresses0[1:])#/stresses0[1:]
	else:
	  diffs = (pystress[1:]-stresses0[1:])/stresses0[1:]
	
	if minmax:
          Fv = np.max([Fv,1000*np.max(np.abs(diffs))])
        else:
          Fv += np.sum(diffs*diffs)/diffs.size
	#Fv += np.sum(diffs*diffs)
	#lstdiffs = np.where(diffs<>np.max(diffs))[0]
        #Fv += np.sum(diffs[lstdiffs]*diffs[lstdiffs])
  
  if printopt:
    printstatem =  "FV: %f  "+", %f"*X.size
    print printstatem%tuple(np.r_[Fv,X])
    
  return Fv
  
 

 
X0 = np.array([V[Xnames[xnr]] for xnr in range(Xnames.__len__())])
SF = 10**(np.log10(X0+10))
X0 /= SF

if doplots:
  p.close(1)
  p.figure(1)
  dtimemax = dtimemax0
  for cnr in range(nrdata):# [1]:
    cnc = cnames[cnr]
    curve = data[cnc]
    times0 = curve[:,0]
    strains0 = curve[:,1]
    stresses0 = curve[:,2]
    temps0 = curve[:,3]
    if cutfigs:
      mxpt = np.where(stresses0==np.max(stresses0))[0][0]
      strains0 = strains0[:mxpt]
      stresses0 = stresses0[:mxpt]
    pystress = getPYstress(strains0,times0,temps0,**V)
    p.plot(strains0,stresses0,'k.',ms=2)#ccolor[cnr]+'.')
    p.plot(strains0,pystress,'k-')#ccolor[cnr])
 
  p.xlabel(r'$\varepsilon$',fontsize=16.)    
  p.ylabel(r'$\sigma$',fontsize=16.)
  p.title('Optimised Parameter Fit')
  p.legend(['Data','Model'],loc='upper left')#'lower right')
  
  if usetests:
    p.close(2)
    p.figure(2)
    dtimemax = dtimemaxtest
    for cnr in range(nrtests):# [1]:
      #p.close(2+cnr)
      #p.figure(2+cnr)
      cnc = testnames[cnr]
      curve = testdata[cnc]
      times0 = curve[:,0]
      strains0 = curve[:,1]
      stresses0 = curve[:,2]
      temps0 = curve[:,3]
      pystress = getPYstress(strains0,times0,temps0,**V)
      p.plot(strains0,stresses0,'k.',ms=2)#ccolor[cnr]+'.')
      p.plot(strains0,pystress,'k-')#ccolor[cnr])
 
      p.xlabel(r'$\varepsilon$',fontsize=16.)    
      p.ylabel(r'$\sigma$',fontsize=16.)
      p.title('Optimised Parameter Fit')
      p.legend(['Data','Model'],loc='upper left')#lower right')
    


time.sleep(2)