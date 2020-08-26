# a modification of code in "git/colr_dep_ESD/esd_plots_obs+theo.py"

import matplotlib.pyplot as plt 
from glob import glob
import numpy as np
from astropy.io import fits
import json
from subprocess import call

def ratioplot_esd(rbin, pofz, method, single):
    #open caluculated esd from aum.
    #fixed dir tree. try not to alter it.---26Aug2020
    base0 = glob(f"/home/navin/git/colr_dep_ESD/{pofz}/run{run}*/signal_dr72safe_bin*")[0]
    with open(f'{base0}/esd_{method}_magbinned_colrdep.json','r') as f: 
        dic = json.load(f)
    magbin = np.array(list(dic.keys()))
    binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
    rp = dic['rp']

    #fixed dir tree. try not to alter it.---26Aug2020
    #open config file
    with open(glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*/*config*")[0], "r") as yamlfile:
        config = yaml.load(yamlfile)
 
    #base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/signal_dr72safe' #firt run --> rbin=15
    base = glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*/signal*")[0]+'/signal_dr72safe'
    sample_base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/col_dep_filter_fits/post_catalog.dr72safe'
    identifier = [7,8,9,10]
    colour = ['red','blue']
    for tag in identifier:

        if single==False:
            #for both colr on the same plot 
            fig,axes = plt.subplots(1,1)
        for colr in colour:
            # get info of samples to be plotted from the red/blue filtered fits files
            hdul = fits.open('%s%d.%s.fits'%(sample_base,tag,colr))
            hdr = hdul[0].header
            leg1 = "abs_mag_bin:(%s,%s)"%(hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3])
            leg2 = "z-bin: (%s)" %(hdr['z_bin'].split('-')[0][:4]+'-'+(hdr['z_bin'].split('-')[1][:4]))
            leg3 = '%s galaxy sample'%colr

            # get data from weaklensing pipeline output to plot
            deltasigma,r,errDeltaSigma = np.loadtxt('%s%d_%s.dat'%(base,tag,colr),usecols=(5,7,12),unpack=True) 

            #to apply filter to rp to get only those radii which are valid points in observed signal files of each color.
            rp_ = np.array(rp)

            #debugging step
            print(f"rp_.size={rp_.size},r.size={r.size}: should be like ({rbin},<={rbin}).")
            print(f"{r==rp_}: if True apply same filter on 'rp' values as on 'r'.")
            print(f"{leg1}\nNo of nan values in \"{colr}\" esd/deltasigma,r,errDeltaSigma: {np.sum(np.isnan(deltasigma))},{np.sum(np.isnan(r))},{np.sum(np.isnan(errDeltaSigma))}")
            print(f"No of negative points in \"{colr}\" esd/deltasigma,errDeltaSigma: {np.sum(deltasigma<0)},{np.sum(errDeltaSigma<0)}" ) #\n\n

            notnan = ~np.isnan(deltasigma) & ~np.isnan(r) & ~np.isnan(errDeltaSigma)
            positives = ~(deltasigma<0) & ~(errDeltaSigma<0)
            deltasigma=deltasigma[notnan & positives]
            r = r[notnan & positives]
            errDeltaSigma = errDeltaSigma[notnan & positives]

            #get calculated(theoretical) ESD from AUM
            if colr=='red':
                esd = dic[magbin[1:][binmax==float(hdr['absmmax'])][0]][0]['red']
            else: 
                esd = dic[magbin[1:][binmax==float(hdr['absmmax'])][0]][1]['blue']

            #to get the ratio plot
            rp_ = rp_[notnan & positives] #get same points of radii at which 'deltasigma' is available. 
            esd = np.array(esd)[notnan & positives]
            #debugging step
            print(f"predicted esd.size={esd.size}\tobserved deltasigma.size{deltasigma.size}\tr(weaklens output files) size={r.size}\trp(ESD_aum)={rp_.size}\n\n")
            ratio = deltasigma/esd #obervational/theoretical

            if single==True:
                # plotting one colr on one plot
                fig,axes = plt.subplots(1,1)
                axes.plot(rp_, ratio, label=r"$\frac{Observed(NYU-VAGC\ cat.+Pipeline)}{Prediction(Niladri\ et\ al.+AUM)}$",c=colr) #marker='d', markerfacecolor='white',
                axes.errorbar([],[],markerfacecolor='None',ls='',label=f"{leg1}\n{leg2}")
                axes.set_xscale('log')
                #axes.set_yscale('log')    
                axes.legend(loc='best', frameon=True)
                plt.title(f"{config['random']} rotation and {config['pofz']} for the HSC sources")
                # setting common x and y axes lables.
                fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
                #fig.text(0.04, 0.5, r'$\frac{Predicted(\Delta\Sigma)}{Observed(\Delta\Sigma)}$', va='center', rotation='vertical')#, fontsize=16) #units_deltasigma=[hM$_\odot$/pc$^2]
              
                # for common title to all the subplots
                plt.savefig((base0 + f"/{method}_{colr}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
            
            if single==False:
                axes.plot(rp_, ratio,c=colr) #marker='d', markerfacecolor='white'
        if single==False:
            plt.title(f"{config['random']} rotation and {config['pofz']} for the HSC sources")
            axes.errorbar([],[],markerfacecolor='None',ls='',label=f"{leg1}\n{leg2}")
            axes.errorbar([],[],markerfacecolor='None',ls='',label=r"$\frac{Observed(NYU-VAGC\ cat.+Pipeline)}{Prediction(Niladri\ et\ al.+AUM)}$")
            axes.set_xscale('log')
            #axes.set_yscale('log')    
            axes.legend(loc='best', frameon=True)
            # setting common x and y axes lables.
            fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
            #fig.text(0.04, 0.5, r'$\frac{Predicted(\Delta\Sigma)}{Observed(\Delta\Sigma)}$', va='center', rotation='vertical')#, fontsize=16)

            ## for common title to all the subplots
            #plt.suptitle(r'samples used for weaklensing single calculation (NYU_VAGC,Zehavi,Niladri)')
            #plt.suptitle('weak lensing signal around NYU-VAGC galaxy sample dr72safe%d in HSC field of sources.\n\n%s, %s'%(tag,leg1,leg2))#, fontsize=15)

            plt.savefig(base0 + f"/{method}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")

if __name__=="__main__":
    #rbin = 5 #same rbin supplied in config file of weaklens_pipeline---used in old version

    method = ['bestfit','fittingFunc']
    pofz = "fullpofz" #"noFullpofz" #fullpofz (give one of two types.)
    #for run in [1,2,3,4,5]: #for noFullpofz
    for run in [0,1,2,3,4,5]: #for fullpofz
        base0 = glob(f"/home/navin/git/colr_dep_ESD/{pofz}/run{run}*/signal_dr72safe_bin*")[0]
        for ii in method:

            for jj,tt in zip((True,False),("single","overplot")):
                ratioplot_esd(run, pofz, method=ii,single=jj) 
                cmd="mkdir -p " + base0 + f"/{tt}.{ii}.ratio"
                call(cmd,shell=True)
                call(f"mv {base0}/*png {base0}/{tt}.{ii}.ratio",shell=True`) 
                call(f"tar -czvf {base0}/{tt}.{ii}.ratio.tar.gz {base0}/{tt}.{ii}.ratio",shell=True) 
        call(f"mkdir -p {base0}/esd_ratio_plots",shell=True)
        call(f"mv -f {base0}/ratio* {base0}/ratio.tar.gz* {base0}/esd_ratio_plots",shell=True)

"""
understand units of deltasigma theoretically as well as the way it gets caclulated in ESd and weaklens pipeline(observational).

removed 'nan' and 'negative' 'deltasigma' points from 'deltasigma' of obserevd weaklensing pipeline outputs. 
correspondingly the 'esd' points from ESD(using AUM and HOD) with same radii are not plotted. 
"""
