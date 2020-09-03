# an extension of code in "Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/plotSignal1.py"  

# variables to play with: base, rbin/ run, pofz 
# the location of json file for esd signal is fixed now...only need to give it 'run' and 'pofz' to access the right file.
# name of config file should always end with "config" ------> NOT USING ANYMORE insdead...I access the ifrandom and pofz from the directory structure itself..which is anyways coming from the actual config file

import matplotlib.pyplot as plt 
from glob import glob
import numpy as np
from astropy.io import fits
import json
from subprocess import call
import yaml

def plot_esd(pofz,rbin, method, single):
    #fixed dir tree. try not to alter it.---26Aug2020, 03Sept2020
    #base0 = glob(f"/home/navin/git/colr_dep_ESD/{pofz}/run{run}*/signal_dr72safe_bin*")[0]
    #base0 = glob(f"/home/navin/git/colr_dep_ESD/{pofz}/run*{run}*{pofz}")[0] #03Sept2020
    base0 = glob(f"/home/navin/git/colr_dep_ESD/Theoretical_ESD_signal")[0] #03Sept2020
    with open(f'{base0}/esd_{method}_magbinned_colrdep_rbin{rbin}.json','r') as f: 
        dic = json.load(f)
    magbin = np.array(list(dic.keys()))
    binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
    rp = dic['rp']
   
    #fixed dir tree. try not to alter it.---26Aug2020
    #open config file
    #no need to store config file for every run. I'll try to extract the required info from the tree structure itself..like below
    #with open(glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*/*config*")[0], "r") as yamlfile:
    #    config = yaml.load(yamlfile)
    #outputdir = glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*")[0] #03Sept2020
    outputdir = glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/rbin{rbin}/run{run}*")[0] #03Sept2020
    ifrandom, pofz = outputdir.split('_')[-3], outputdir.split('_')[-1] 
    
    #base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/signal_dr72safe'
    #base = glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*/signal*")[0]+'/signal_dr72safe' #26Aug2020
    base = glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*")[0]+'/signal_dr72safe' #03Sept2020
    #Specific to downloaded nyu_data..No need to change.
    sample_base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/col_dep_filter_fits/post_catalog.dr72safe'
    identifier = [7,8,9,10]
    colour = ['red','blue']
    for tag in identifier:
        #fig,axes = plt.subplots(nrows=2,ncols=2)
        #for ii,colr in zip(range(len(axes)),colour): #accessing one row at a time.
        
        if single==False:
            #for both colr on the same plot 
            fig,axes = plt.subplots(1,1)
        for colr in colour:
            # get data from weaklensing pipeline output to plot
            if ifrandom == 'random':
            deltasigma,r,errDeltaSigma = np.loadtxt('%s%d_%s.dat'%(base,tag,colr),usecols=(5,7,12),unpack=True)
            notnan = ~np.isnan(deltasigma) & ~np.isnan(r) & ~np.isnan(errDeltaSigma)
            deltasigma=deltasigma[notnan]
            r = r[notnan]
            errDeltaSigma = errDeltaSigma[notnan]
            # get info of samples to be plotted from the red/blue filtered fits files
            hdul = fits.open('%s%d.%s.fits'%(sample_base,tag,colr))
            hdr = hdul[0].header
            leg1 = "abs_mag_bin:(%s,%s)"%(hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3])
            leg2 = "z-bin: %s" %hdr['z_bin'].split('-')[0][:4]+'-'+(hdr['z_bin'].split('-')[1][:4])
            leg3 = '%s galaxy sample'%colr
            print("Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):")
            print(f"{leg1}\nNo of nan values in \"{colr}\" esd/deltasigma,r,errDeltaSigma: {np.sum(np.isnan(deltasigma))},{np.sum(np.isnan(r))},{np.sum(np.isnan(errDeltaSigma))}")
            print(f"No of negative points in\"{colr}\" esd/deltasigma,errDeltaSigma: {np.sum(deltasigma<0)},{np.sum(errDeltaSigma<0)}\n\n" )

            #get calculated(theoretical) ESD from AUM
            if colr=='red':
                esd = dic[magbin[1:][binmax==float(hdr['absmmax'])][0]][0]['red']
            else: 
                esd = dic[magbin[1:][binmax==float(hdr['absmmax'])][0]][1]['blue']
            if single==True:
                # plotting one colr on one plot
                fig,axes = plt.subplots(1,1)
                axes.plot(rp, esd, label=f"Prediction:(model HOD)Niladri et al.",c=colr) #marker='d', markerfacecolor='white',
                axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC_data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
                axes.errorbar([],[],markerfacecolor='None',ls='',label=f"{leg1}\n{leg2}")
                #axes.scatter([],[],facecolors='None',label="=======")
                axes.set_ylim(0.01,3000)
                axes.set_xscale('log')
                axes.set_yscale('log')    
                axes.legend(loc='best', frameon=True)
                plt.title(f"{ifrandom} rotation and {pofz} for the HSC sources")
                #plt.title(f"{method} params used to get HODs.(Niladri et al.)")
                # setting common x and y axes lables.
                fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
                fig.text(0.04, 0.5, '$\Delta\Sigma$ [hM$_\odot$/pc$^2]$', va='center', rotation='vertical')#, fontsize=16)
                # for common title to all the subplots
                plt.savefig(base0 + f"/{method}_{colr}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
            else:
                # plotting red and blue on the same plot 
                if colr=='blue':
                    axes.plot(rp, esd, label=f"Prediction:(HOD model)Niladri et al.",c=colr) #marker='d', markerfacecolor='white'
                    axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC_data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
                else:
                    axes.plot(rp, esd, c=colr) #, marker='d',markerfacecolor='white'
                    axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', fmt='o', linewidth=1.5, capsize=5, capthick=2)
        if single==False:
            plt.title(f"{ifrandom} rotation and {pofz} for the HSC sources")
            axes.errorbar([],[],markerfacecolor='None',ls='',label=f"{leg1}\n{leg2}")
            axes.set_ylim(0.01,3000)
            axes.set_xscale('log')
            axes.set_yscale('log')    
            axes.legend(loc='best', frameon=True)
            # setting common x and y axes lables.
            fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
            fig.text(0.04, 0.5, '$\Delta\Sigma$ [hM$_\odot$/pc$^2]$', va='center', rotation='vertical')#, fontsize=16)
            #plt.title(f"{method} params used to get HODs.(Niladri et al.)")
            # for common title to all the subplots
            #plt.suptitle(r'samples used for weaklensing signal calculation (NYU_VAGC,Zehavi,Niladri)')
            #plt.suptitle('weak lensing signal around NYU-VAGC galaxy sample dr72safe%d in HSC field of sources.\n\n%s, %s'%(tag,leg1,leg2))#, fontsize=15)
            #plt.show()
            plt.savefig(base0 + f"/{method}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
   
if __name__=="__main__":
    method = ['bestfit','fittingFunc']
    #rbin = 10 #5 #15  #change it as per the config file for weaklens_pipeline---older version.

    # run and pofz decide the location of files to be plotted and stored---26Aug2020
    #run = 2
    pofz = "nofullpofz" #"nofullpofz" #fullpofz (give one of two types.)
    #for run in [0,1,2,3,4,5]: #for nofullpofz
    for run in [0,1,2,3,4,5]: #for fullpofz
        #base0 = glob(f"/home/navin/git/colr_dep_ESD/{pofz}/run{run}*/signal_dr72safe_bin*")[0]
        base0 = glob(f"/home/navin/git/colr_dep_ESD/{pofz}/run{run}*")[0]
        for ii in method:
            for jj,tt in zip((True,False),("single","overplot")):
                plot_esd(run, pofz, method=ii, single=jj) 
                cmd = "mkdir -p " + base0 + f"/{tt}.{ii}"
                call(cmd,shell=True)
                call(f"mv {base0}/*png {base0}/{tt}.{ii}",shell=True)
                call(f"tar -czvf {base0}/{tt}.{ii}.tar.gz {base0}/{tt}.{ii}",shell=True)
        call(f"mkdir -p {base0}/esd_plots",shell=True)
        call(f"mv -f {base0}/single* {base0}/overplot* {base0}/esd_plots",shell=True)



"""

At some stage I check in the files outputted from weaklens pipeline:
1)There are some negative values in the esd/deltasigma in the observed signal(which was calulated from weaklens pipeline.)
Do you understan why would esd be negative theoretically?
2) Also there are nan values in deltasigma and errDeltaSigma.
-------> nan values are filtered out but negative values are not yet removed...
------->-ve values doesn't cause any problem in plotting since we are not taking log()of these numbers...only scale is in log.
------->But in ratio plots these -ve points are excluded.

Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):
Notice: nan values are absent only because the prints stements are taken after filter nan entries.
abs_mag_bin:(-23.0,-22.0)
No of nan values in "red" esd/deltasigma,r,errDeltaSigma: 0,0,0
No of negative points in"red" esd/deltasigma,errDeltaSigma: 0,0


Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):
abs_mag_bin:(-23.0,-22.0)
No of nan values in "blue" esd/deltasigma,r,errDeltaSigma: 0,0,0
No of negative points in"blue" esd/deltasigma,errDeltaSigma: 0,0


Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):
abs_mag_bin:(-22.0,-21.0)
No of nan values in "red" esd/deltasigma,r,errDeltaSigma: 0,0,0
No of negative points in"red" esd/deltasigma,errDeltaSigma: 0,0


Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):
abs_mag_bin:(-22.0,-21.0)
No of nan values in "blue" esd/deltasigma,r,errDeltaSigma: 0,0,0
No of negative points in"blue" esd/deltasigma,errDeltaSigma: 0,0


Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):
abs_mag_bin:(-21.0,-20.0)
No of nan values in "red" esd/deltasigma,r,errDeltaSigma: 0,0,0
No of negative points in"red" esd/deltasigma,errDeltaSigma: 0,0


Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):
abs_mag_bin:(-21.0,-20.0)
No of nan values in "blue" esd/deltasigma,r,errDeltaSigma: 0,0,0
No of negative points in"blue" esd/deltasigma,errDeltaSigma: 7,0


Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):
abs_mag_bin:(-20.0,-19.0)
No of nan values in "red" esd/deltasigma,r,errDeltaSigma: 0,0,0
No of negative points in"red" esd/deltasigma,errDeltaSigma: 3,0


Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):
abs_mag_bin:(-20.0,-19.0)
No of nan values in "blue" esd/deltasigma,r,errDeltaSigma: 0,0,0
No of negative points in"blue" esd/deltasigma,errDeltaSigma: 6,0
"""
