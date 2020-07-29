# Use hod.ESD() from hod.py to get the lesing signal.
#
#For documentation on 'aum' and inside modules, either visit http://surhudm.github.io/aum/code.html or try doing dir(),help(), and ? and ?? in front of every python object defined inside these modules.

## pre-requisites to use this function:
   # there are defined functions in this module, to get c-array from a numpy array  and vice-versa. These conversions of objects types are needed inside ESD().
   # __init__ initializes hod and cosmology object. we'll have to use it to initialize the hod and cosmology class.
   # You might want to use __swig_destroy__() when you are doing some calculation which includes more than one cosmology.
   # Some functions which might be of interest:
      # getOmb, getOmk, geth,gets8, hod_free, hod_renew, ncen, ncenz, nsat, nsatz, set_cen_offset_paramsm, set_cfactor, set_inc_params, sethalo_exc
#__________________________________

"""Plan"""
# set TINK==2 in hod.h (delete the ".../aum/build" dir and then run "python setup.py install --prefix=`pwd`/install" twice in the ".../aum/" dir. This create the ".../aum/build" dir again.)
# import red/blue HODs which match Niladri's best fit results of Global Analysis in numpy arrays.
# convert binned logM and Nsat/Ncen (red/blue) numpy arrays into c-arrays.
# supply these arrays appropriately to set the gsl_interpolation in hod.cpp
# set up other input parameters to call wp_ESD func from hod.py
# save and plot the wp_ESD signals for red/blue HODs...while using (best -fit) and (fitting func) separately...and also compare these two....Discuss with Surhud and Group!!  

"""Manual supply:""" 
# sampled_hod_loc
# cosmo() initialization parameters
#__________________________________

import sys
sys.path.append("/home/navin/aum/install/lib/python3.7/site-packages/")
import hod as h
import cosmology as cc
import numpy as np
from subprocess import call
from glob import glob

# getting c-array from a np array
def getdblarr(r):
    temp=h.doubleArray(r.size)
    for i in range(r.size):
        temp[i]=r[i]
    return temp
# getting np array from a c-array
def getnparr(r,n):
    temp=np.zeros(n)
    for i in range(n):
        temp[i]=r[i]
    return temp



# initialize cosmology class object
# the following cosmo parameters are taken from Niladri's paper
# what are w0, wa, ximax? 
# the hod params are supplied just to initialize the hod object.
# these params won't be used anyways since to get ncen and nsat, interpolation will be used.(TINK==2) 
p = h.cosmo()
q = h.hodpars()
p.Om0 = 0.276
p.w0 = -1
p.wa = 0
p.Omk = 0.0
p.hval = 0.7
p.Omb = 0.045
p.th = 2.726
p.s8 = 0.811
p.nspec = 0.961
p.ximax = np.log10(8.0)
p.cfac = 1.0
q.Mmin = 13.0 
q.siglogM = 0.5 
q.Msat = 14.0 
q.alpsat = 1.0 
q.Mcut = 13.5 
q.csbycdm = 1.0 
q.fac = 1.0
 
a = h.hod(p, q) 

#Mr_max = np.arange(-19,-23,-1)
## number of bins np.size(Mr_max)
#bins = np.array(['%d,%d'%(l-1,l) for l in Mr_max ])

# first write code to get only one file then we'll import other files in a loop based on magbin and colr.
# for now I go  ahead with red HODs.

#set up sampled hod locations
galtype=['cen','sat']
colr=["red","blue"] #use glob to get only those files containing 'red'...
sampled_hod_loc = "/home/navin/git/hod_red_blue/bestfit_binned_hods/"
cen_hod_loc, sat_hod_loc = [sampled_hod_loc+ x for x in galtype]

#store file_names based on colr-galtype in increasing brightness oreder 
temp0 = [glob(cen_hod_loc+f"/*{x}*") for x in colr] 
temp1 = [glob(sat_hod_loc+f"/*{x}*") for x in colr]
[x.sort() for x in temp0]
[x.sort() for x in temp1]

# Note: the negative values in binned_colr dep hods---(Ncen_red,Ncen_blue,Nsat_red,Nsat_blue) 
# a bug in model??
for ii,col in enumerate(colr):
    for jj,colr_pair in enumerate(zip(temp0[ii],temp1[ii])):
        print(ii,col,jj,colr_pair)
        #cen hod
        logM, hod0 = np.loadtxt(colr_pair[0],dtype={'names':("logM","hod",), 'formats': ('float','float',)},comments="#", unpack=True)
        #sat hod
        _ , hod1 = np.loadtxt(colr_pair[1],dtype={'names':("logM","hod",), 'formats': ('float','float',)},comments="#", unpack=True)
        # initialize spline, TINK==2
        a.init_Nc_spl(getdblarr(logM[hod0>0]), getdblarr(np.log10(hod0[hod0>0])), hod0[hod0>0].size)
        a.init_Ns_spl(getdblarr(logM[hod1>0]), getdblarr(np.log10(hod1[hod1>0])), hod1[hod1>1].size)
        #debug step
        print(f"{col}_hod,for mass in {np.arange(11.0,16.0,1.0)}\nncen={list( map(a.ncen,np.arange(11.0,16.0,1.0)))}\n nsat={list(map(a.nsat,np.arange(11.0,16.0,1.0)))}")
