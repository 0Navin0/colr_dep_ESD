# an extension of code in "Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/plotSignal1.py"  

# variables to play with: base, rbin 

import matplotlib.pyplot as plt 
from glob import glob
import numpy as np
from astropy.io import fits
import json
from subprocess import call

def plot_esd(rbin, method, single):
    with open(f'./bin{rbin}/esd_{method}_magbinned_colrdep.json','r') as f: 
        dic = json.load(f)
    magbin = np.array(list(dic.keys()))
    binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
    rp = dic['rp']
    
    #base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/signal_dr72safe'
    base = f'/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/signal_dr72safe_bin{rbin}/signal_dr72safe'
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
                #plt.title(f"{method} params used to get HODs.(Niladri et al.)")
                # setting common x and y axes lables.
                fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
                fig.text(0.04, 0.5, '$\Delta\Sigma$ [hM$_\odot$/pc$^2]$', va='center', rotation='vertical')#, fontsize=16)
                # for common title to all the subplots
                plt.savefig(f"{method}_{colr}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
            else:
                # plotting red and blue on the same plot 
                if colr=='blue':
                    axes.plot(rp, esd, label=f"Prediction:(HOD model)Niladri et al.",c=colr) #marker='d', markerfacecolor='white'
                    axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC_data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
                else:
                    axes.plot(rp, esd, c=colr) #, marker='d',markerfacecolor='white'
                    axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', fmt='o', linewidth=1.5, capsize=5, capthick=2)
        if single==False:
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
            plt.savefig(f"{method}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
   
if __name__=="__main__":
    method = ['bestfit','fittingFunc']
    rbin = 10 #5 #15  #change it as per the config file for weaklens_pipeline
    for ii in method:
        for jj,tt in zip((True,False),("single","overplot")):
            plot_esd(rbin, method=ii, single=jj) 
            cmd=f"mkdir ./bin{rbin}/{tt}.{ii}"
            call(cmd,shell=True)
            call(f"mv *png ./bin{rbin}/{tt}.{ii}",shell=True)
            call(f"tar -czvf ./bin{rbin}/{tt}.{ii}.tar.gz ./bin{rbin}/{tt}.{ii}",shell=True)
    call(f"mkdir -p ./bin{rbin}/esd_plots",shell=True)
    call(f"mv ./bin{rbin}/single* ./bin{rbin}/overplot* ./bin{rbin}/esd_plots",shell=True)



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
