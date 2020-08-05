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
# set TINK==2 in hod.h (delete the ".../aum/build" dir and then run "python setup.py install --prefix=`pwd`/install" twice in the ".../aum/" dir. This creates the ".../aum/build" dir again.)
# import red/blue HODs which match Niladri's best fit results of Global Analysis in numpy arrays. (/home/navin/git/hod_red_blue/bestfit_binned_hods/colr)
# convert binned logM and Nsat/Ncen (red/blue) numpy arrays into c-arrays.
# supply these arrays appropriately to set the gsl_interpolation in hod.cpp
# set up other input parameters to call wp_ESD func from hod.py
# save and plot the wp_ESD signals for red/blue HODs...while using (best-fit) and (fitting func) separately...and also compare these two....Discuss with Surhud and Group!!  

"""Manual supply:""" 
# sampled_hod_loc
# cosmo() initialization parameters
#__________________________________

#accessing aum
import sys
sys.path.append("/home/navin/aum/install/lib/python3.7/site-packages/")
import hod as h
import cosmology as cc
#dealing with files in local system 
from subprocess import call
from glob import glob
#manipulation of data
import numpy as np
import pandas as pd
#managing and saving data
from collections import defaultdict 
import json

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

def initializeHOD():   
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
     
    return h.hod(p, q) 


def esd_json(method='bestfit'):
    
    #set up sampled hod locations
    galtype=['cen','sat']
    colr=["red","blue"] #use glob to get only those files containing 'red'...
    sampled_hod_loc = "/home/navin/git/hod_red_blue/%s_binned_hods/"% method
    cen_hod_loc, sat_hod_loc = [sampled_hod_loc+ x for x in galtype]
    
    #store HOD file_names based on colr-galtype in increasing brightness order 
    temp0 = [glob(cen_hod_loc+f"/*{x}*") for x in colr] 
    temp1 = [glob(sat_hod_loc+f"/*{x}*") for x in colr]
    [x.sort() for x in temp0]
    [x.sort() for x in temp1]
    #get the number of magbins and magbins to match the correct redshift to each sample.
    mag_order=[temp0[0][jj].split('_')[-1].split('.csv')[0] for jj in range(len(temp0[0]))]
    
    #get mean redshift from the samples(pr.magbin.dat) in the order same as mag_order
    df = pd.read_csv("/home/navin/git/hod_red_blue/pr.magbin.dat", delim_whitespace=True)
    #z_order1=[(f"{df.mag1[ii]}"+"-"+f"{df.mag2[ii]}",df.z[ii]) for ii in range(df.mag2.size)]
    z_order=[f"{df.mag1[ii]}"+"-"+f"{df.mag2[ii]}" for ii in range(df.mag2.size) ]
    df['z_order']=z_order
    
    #define proj-radii and esdbins for ESD caclucation from aum
    #the radii used in observed esd by weaklens pipeline
    base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/signal_dr72safe'
    rp = np.loadtxt('%s7_red.dat'%base,usecols=(7,),unpack=True) #same for all magbins and colr
    #rp =np.linspace(0.02,4,15) 
    esdbins = rp.size
    
    #store all the weak lensing signal data.
    res = defaultdict(list)
    #ndarray objects can't be serialised to json object..so convert wls to list(rp)
    res['rp']=list(rp)    
    #get aum ready
    a = initializeHOD()

    # Note: the negative values in binned_colr dep hods---(Ncen_red,Ncen_blue,Nsat_red,Nsat_blue) 
    # a drawback of model?? --> Yes.These negative values are unphysical.
    for ii,col in enumerate(colr):
        for jj,colr_pair in enumerate(zip(temp0[ii],temp1[ii])):
            ## get redshift of the galaxy sample(NYU catalog safe7 sample) in magbin-mag_order
            #z = next((z_order1[ii][1] for ii in range(len(z_order1)) if mag_order[jj]==z_order1[ii][0])) #another way with z_order1
            z = df.z.values[df['z_order']==mag_order[jj]][0]
            print(ii,col,jj,colr_pair,mag_order[jj],z)
    
            ## cen hod
            logM, hod0 = np.loadtxt(colr_pair[0],dtype={'names':("logM","hod",), 'formats': ('float','float',)},comments="#", unpack=True)
            ## sat hod
            _ , hod1 = np.loadtxt(colr_pair[1],dtype={'names':("logM","hod",), 'formats': ('float','float',)},comments="#", unpack=True)
            print(f"cen:{hod0[hod0>0].size}, sat:{hod1[hod1>0].size}")
            hod0[hod0<=0]=1e-50 
            hod1[hod1<=0]=1e-50
            print(f"cen:{hod0.size}, sat:{hod1.size}")
            #print(f"cen:{hod0},sat:{hod1}")
 
            ## initialize spline, TINK==2
            a.init_Nc_spl(getdblarr(logM), getdblarr(np.log10(hod0)), hod0.size)
            a.init_Ns_spl(getdblarr(logM), getdblarr(np.log10(hod1)), hod1.size)
            ## debug step
            print(f"{col}_hod,for mass in {np.arange(11.0,16.0,1.0)}\nncen={list( map(a.ncen,np.arange(11.0,16.0,1.0)))}\n nsat={list(map(a.nsat,np.arange(11.0,16.0,1.0)))}")
            #print(f"{esdbins}, projected radii={rp}")
            esdrp = getdblarr(rp)
            esd = getdblarr(np.zeros(esdbins))
            a.ESD(z,esdbins,esdrp,esd,esdbins+4)
            wls = getnparr(esd,esdbins)
            print(col,mag_order[jj],wls)
            ##ndarray objects can't be serialised to json object..so convert wls to list(wls)
            res[mag_order[jj]].append({col:list(wls)}) 
    
    with open(f"esd_{method}_magbinned_colrdep.json", "w") as f:
         json.dump(res, f)

if __name__=="__main__":
    esd_json(method="bestfit")     
    esd_json(method="fittingFunc")     







"""
	## to decerialise the json file into a dict:
        with open('esd_binned_colrdep.json','r') as f: 
            dic = json.load(f) 
"""

"""
	accessing red and blue EDDs from dictionary res/dic
	dic['-19.0--20.0'][0]['red'] ,dic['-20.0--21.0'][0]['red']
	dic['-19.0--20.0'][1]['blue'],dic['-20.0--21.0'][1]['blue']
"""

"""
        #confirm the z values supplied for each magbin and col - done (approximated values are entered for ESD calculation.)

#red_mag_bin	z_min to z_max	 	N_gal		N_red		N_blue		safe_sample from NYU(SDSS)		average 'z'(red sample)		average 'z'(blue sample)	average 'z'(red+blue) 
#(-23,-22)	0.1000000-0.2500000	10857		8926		1931  		post_catalog.dr72safe7			0.19119226454406368 		0.18408367486891508 		0.18763797
#(-22,-21)	0.0663000-0.1588300	73628		46021		27607		post_catalog.dr72safe8			0.12267274637529606		0.12368363747192741		0.12317819	
#(-21,-20)*	0.0420000-0.0641700	17825		9807		8018		post_catalog.dr72safe9			0.0545910274847526 		0.05483306938051494 		0.05471205 
#(-20,-19)	0.0268300-0.0641700	44264		19018		25246		post_catalog.dr72safe10			0.04946038991021664		0.05003913348123218		0.04974976 

        #confirm the rp and esdbins supplied. - done.
        #next we save the ESD for col and magbin in some file. - done
        #setting hod_value = 1e-50 may cause ringing effects. Pay attention to this issue. ---- ?
"""
