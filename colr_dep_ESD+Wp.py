# Use hod.wp_ESD() from hod.py to get the lesing signal.(projected Excess Surface Density)
#
#For documentation on 'aum' and inside modules, either visit http://surhudm.github.io/aum/code.html or try doing dir(),help(), and ? and ?? in front of every python object defined inside these modules.

## pre-requisites to use this function:
   # there are defined functions in this module, to get c-array from a numpy array  and vice-versa. These conversions of objects types are needed inside ESD().
   # we need to initialize the cosmology and hod class.
   # You might want to use __swig_destroy__() when you are doing some calculation which includes more than one cosmology.
   # Some functions which might be of interest:
      # getOmb, getOmk, geth,gets8, hod_free, hod_renew, ncen, ncenz, nsat, nsatz, set_cen_offset_paramsm, set_cfactor, set_inc_params, sethalo_exc
#__________________________________

"""Plan"""
# set TINK==2 in hod.h (delete the ".../aum/build" dir and then run "python setup.py install --prefix=`pwd`/install" twice in the ".../aum/" dir. This creates the ".../aum/build" dir again.)
# import red/blue HODs which match Niladri's best fit results of Global Analysis in numpy arrays. (/home/navin/git/hod_red_blue/bestfit_binned_hods/cen(or sat))
# convert binned logM and Nsat/Ncen (red/blue) numpy arrays into c-arrays.
# supply these arrays appropriately to set the gsl_interpolation in hod.cpp
# set up other input parameters to call wp_ESD func from hod.py
# save and plot the wp_ESD signals for red/blue HODs...while using (best-fit) and (fitting func) separately...and also compare these two....Discuss with Surhud and Group!!  
# similaryly calculate Wp_rp (clustering signal) in each bin for each color(red/blue).

"""Manual supply:""" 
# sampled_hod_loc
# cosmo() initialization parameters
# path to 'base' variable.(to get the same radial binning as in the weaklens_pipeline for the observed signals.)
# rbin or run, pofz
# __________________________________

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
#used to open the config file to read some parameters for labelling and naming the files
import yaml

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


def esd_json(rbin, get_wp=False, get_esd=True, method='bestfit', save_interp_hod=True):
    
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
    #get the number of magbins and magbins to match the correct redshift to each colr-galtype sample.
    mag_order=[temp0[0][jj].split('_')[-1].split('.csv')[0] for jj in range(len(temp0[0]))]
    
    #get mean redshift for the samples(pr.magbin.dat) in the order same as mag_order
    df = pd.read_csv("/home/navin/git/hod_red_blue/pr.magbin.dat", delim_whitespace=True)
    #z_order1=[(f"{df.mag1[ii]}"+"-"+f"{df.mag2[ii]}",df.z[ii]) for ii in range(df.mag2.size)]
    z_order=[f"{df.mag1[ii]}"+"-"+f"{df.mag2[ii]}" for ii in range(df.mag2.size) ]
    df['z_order']=z_order
   
    if get_esd: 
        #define proj-radii and esdbins for ESD caclucation from aum
        #the same radii used in observed esd by weaklens pipeline
        #base = '/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/signal_dr72safe_bin{rbin}/signal_dr72safe' #17Aug2020
        #base = f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*/signal_dr72safe*"  #26Aug2020 # don't keep produced signals in different rbins in same "run#" dir.
        #base = f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/*fullpofz/rbin{rbin}/run*"  #03Sept2020
        base = f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/*/run*rbin{rbin}*" 
        base = glob(base)[0]
        #debug step
        print(f"\n\nUsing {base} for projected radii bins.\n\n")
        rp = np.loadtxt('%s/signal_dr72safe7_red.dat'%base,usecols=(7,),unpack=True) #same for all magbins and colr
        #store all the weak lensing signal data.
        res = defaultdict(list)#esd
        #ndarray objects can't be serialised to json object..so convert wls to list(rp)
        res['rp']=list(rp)
        esdbins = rp.size  
    """
        #open config file
        with open('/'.join(base.split('/')[:-1]), "r") as yamlfile:
            config = yaml.load(yamlfile)
    """

    if get_wp:
        res1 = defaultdict(list)#wp
        #ndarray objects can't be serialised to json object..so convert np.array(rp) to list(rp)
        res1['rp']=list(np.logspace(-1,np.log10(30),50)) 
        wpbins = len(res1['rp'])

    #get aum ready
    a = initializeHOD()

    # Note: the negative values in binned_colr dep hods---(Ncen_red,Ncen_blue,Nsat_red,Nsat_blue) 
    # a drawback of model?? --> Yes.These negative values are unphysical.
    for ii,col in enumerate(colr):
        for jj,colr_pair in enumerate(zip(temp0[ii],temp1[ii])):
            ## get redshift of the galaxy sample(NYU catalog safe7 sample) in magbin-mag_order
            #z = next((z_order1[ii][1] for ii in range(len(z_order1)) if mag_order[jj]==z_order1[ii][0])) #another way with z_order1
            z = df.z.values[df['z_order']==mag_order[jj]][0]
            print("\n*-----------*")
            print(ii,col,jj,colr_pair,mag_order[jj],z)
 
            ## cen hod
            logM, hod0 = np.loadtxt(colr_pair[0],dtype={'names':("logM","hod",), 'formats': ('float','float',)},comments="#", unpack=True)
            ## sat hod
            _ , hod1 = np.loadtxt(colr_pair[1],dtype={'names':("logM","hod",), 'formats': ('float','float',)},comments="#", unpack=True)
            print(f"fraction of positive points in modeled HODs:\ncen:{hod0[hod0>0].size}/{hod0.size}, sat:{hod1[hod1>0].size}/{hod1.size}")
            print(f"""Negatives and zeros are to be removed. Zeros give divisionbyzeroError, negatives are unphysical. 
In some mag bins negatives come for low mass halos while for some for high mass halos.
Go and check: /home/navin/git/hod_red_blue/{method}_binned_hods/cen(or sat)""")

            #idx = (hod0>0) & (hod1>0)
            #logM = logM[idx]
            #hod0 = hod0[idx]
            #hod1 = hod1[idx]
           
            #removing the zeros and negatives -> cen and sat have different number of points for interpolation. 
            idx0 = hod0>0
            idx1 = hod1>0
            hod0 = hod0[idx0]
            hod1 = hod1[idx1]
            print(f"cen:{hod0.size}, sat:{hod1.size}, logM:{logM[idx0].size}, {logM[idx1].size}")
            #print(f"cen:{h0},sat:{h1}")  
 
            ## initialize spline, TINK==2
            a.init_Nc_spl(getdblarr(logM[idx0]), getdblarr(np.log10(hod0)), hod0.size)
            a.init_Ns_spl(getdblarr(logM[idx1]), getdblarr(np.log10(hod1)), hod1.size)
            ## debug step
            print(f"logM\tcen: {min(logM[idx0])}-{max(logM[idx0])}\t sat: {min(logM[idx1])}-{max(logM[idx1])}")
            #print("*-----------*")
            print(f"{col}_hod,for mass in {np.arange(11.0,16.0,1.0)}\nncen={list( map(a.ncen,np.arange(11.0,16.0,1.0)))}\n nsat={list(map(a.nsat,np.arange(11.0,16.0,1.0)))}")

            if save_interp_hod:
                print("\n.\n.\n.\nsaving the interpolated HODs for testing purposes with 500 interpolated points.\n")
                np.savetxt(colr_pair[0].split('/')[-1], np.c_[np.arange(10.0,16.0,0.01), np.asarray(list( map(a.ncen,np.arange(10.0,16.0,0.01))))], delimiter=' ', header=f"logM {'hod_'+'_'.join((colr_pair[0].split('/')[-1]).split('_')[0:2])}", comments='#')   
                np.savetxt(colr_pair[1].split('/')[-1], np.c_[np.arange(10.0,16.0,0.01), np.asarray(list( map(a.nsat,np.arange(10.0,16.0,0.01))))], delimiter=' ', header=f"logM {'hod_'+'_'.join((colr_pair[1].split('/')[-1]).split('_')[0:2])}", comments='#')   
            print(f"avmass_cen={a.avmass_cen(z)},avmass_tot={a.avmass_tot(z)}")
            
            if get_esd:
                #print(f"esdbins:{esdbins}, projected radii={rp}")
                esdrp = getdblarr(rp)
                esd = getdblarr(np.zeros(esdbins))
                a.ESD(z,esdbins,esdrp,esd,esdbins+4)
                wls = getnparr(esd,esdbins)
                print(col,mag_order[jj],'wl_signal:',wls)
                ##ndarray objects can't be serialised to json object..so convert wls to list(wls)
                res[mag_order[jj]].append({col:list(wls)})

            if get_wp:
                wp_rp = getdblarr(np.array(res1['rp']))
                wp = getdblarr(np.zeros(wpbins))
                a.Wp(z, wpbins, wp_rp, wp, 100.0)
                cls_sig = getnparr(wp,wpbins)
                res1[mag_order[jj]].append({col:list(cls_sig)})
    if get_esd:
        call("mkdir -p ./Theoretical_ESD_signal", shell=True) #03Sept2020
        #call("mkdir -p ./%s"% ('/'.join(base.split('/')[-3:])), shell=True) #26Aug2020
        #call("mkdir -p ./bin{rbin}/{pofz}",shell=True) #17Aug2020 #older version
        #with open(f"./bin{rbin}/{pofz}/esd_{method}_magbinned_colrdep.json", "w") as f:
        #with open(f"./%s/esd_{method}_magbinned_colrdep.json"%('/'.join(base.split('/')[-3:])), "w") as f:
        with open(f"./Theoretical_ESD_signal/esd_{method}_magbinned_colrdep_rbin%s.json"%rbin, "w") as f:
             json.dump(res, f)
             #print(f"esd file saved in: ./%s/esd_{method}_magbinned_colrdep.json\n\n\n"%('/'.join(base.split('/')[-3:])))
             print(f"esd file saved in: ./Theoretical_ESD_signal/esd_{method}_magbinned_colrdep_rbin%s.json\n\n\n"%rbin)
    if get_wp:
        call(f"mkdir -p ./wp",shell=True)
        with open(f"./wp/Wp_{method}_magbinned_colrdep.json", "w") as f:
             json.dump(res1, f)
             print(f"wp file saved in: ./wp/Wp_{method}_magbinned_colrdep.json\n\n\n")

if __name__=="__main__":
    #In newer version, rbin is not needed. --26Aug2020
    #rbin = 10 #used only in address to store esd, same as rbin in config file of weaklens_pipeline
    
    #17Aug2020
    #esd_json(rbin, pofz, get_wp=True, get_esd=False, method="bestfit")     
    #esd_json(rbin, pofz, get_wp=True, get_esd=False, method="fittingFunc")
    
    #26Aug2020 
    #pofz = "fullpofz"# "noFullpofz" #"fullpofz"# (give one of two types.)
    ##for run in [1,2,3,4,5]: #for noFullpofz
    #for run in [0,1,2,3,4,5]: #for fullpofz
    #    #26Aug2020
    #    esd_json(run, pofz, get_wp=False, get_esd=True, method="bestfit")     
    #    esd_json(run, pofz, get_wp=False, get_esd=True, method="fittingFunc")     
   
    for rbin,interp in zip([3,5,10],[False,False,False]):  
        #03Sept2020
        #giving True to any one bin scheme will do since HOD interpolation uses the same data points for all bin sizes from hod_red_blue dir.  
        esd_json(rbin, get_wp=True, get_esd=False, method="bestfit", save_interp_hod=interp)     
        esd_json(rbin, get_wp=True, get_esd=False, method="fittingFunc", save_interp_hod=interp)     






"""
	## to decerialise the json file into a dict:
        with open('esd_binned_colrdep.json','r') as f: 
            dic = json.load(f) 
"""

"""
	accessing red and blue ESDs from dictionary res/dic
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
		yes it causes wrinkles in the plot..I changed this value to 1e-30, 1e-20...you can notice the difference.
		Surhud says, he would never have set 1e-20, but just remove those because AUM will autometically set hod=0.0 there where(for some logM values.) it's not available for interplotation.
"""
