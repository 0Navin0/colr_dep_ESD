# an extension of code in "Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/plotSignal1.py"  

import matplotlib.pyplot as plt 
from glob import glob
import numpy as np
from astropy.io import fits
import json

with open('./esd_bestfit_magbinned_colrdep.json','r') as f: 
    dic = json.load(f)
magbin = np.array(list(dic.keys()))
binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
rp = dic['rp']

base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/signal_dr72safe'
sample_base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/col_dep_filter_fits/post_catalog.dr72safe'
identifier = [7,8,9,10]
colour = ['red','blue']
for tag in identifier:
    #fig,axes = plt.subplots(nrows=2,ncols=2)
    #for ii,colr in zip(range(len(axes)),colour): #accessing one row at a time.
    #fig,axes = plt.subplots(1,1)
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
        #get calculated(theoretical) ESD from AUM
        if colr=='red':
            esd = dic[magbin[1:][binmax==float(hdr['absmmax'])][0]][0]['red']
        else: 
            esd = dic[magbin[1:][binmax==float(hdr['absmmax'])][0]][1]['blue']
     
        # when plotting one colr on one plot
        fig,axes = plt.subplots(1,1)
        axes.plot(rp, esd, label=f"prediction(by AUM): Niladri et al.",c=colr) #marker='d', markerfacecolor='white',
        axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC_data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
        axes.errorbar([],[],markerfacecolor='None',ls='',label=f"{leg1}\n{leg2}")
        #axes.scatter([],[],facecolors='None',label="=======")
        axes.set_xscale('log')
        axes.set_yscale('log')    
        axes.legend(loc='best', frameon=True)
        # setting common x and y axes lables.
        fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
        fig.text(0.04, 0.5, '$\Delta\Sigma$ [hM$_\odot$/pc$^2]$', va='center', rotation='vertical')#, fontsize=16)
        # for common title to all the subplots
        plt.savefig(f"{colr}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")

    #    # plotting red and blue on the same plot (comment above block and uncomment the axes outside this loop.) 
    #    if colr=='blue':
    #        axes.plot(rp, esd, label=f"prediction(by AUM): Niladri et al.",c=colr) #marker='d', markerfacecolor='white'
    #        axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC_data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
    #    else:
    #        axes.plot(rp, esd, c=colr) #, marker='d',markerfacecolor='white'
    #        axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', fmt='o', linewidth=1.5, capsize=5, capthick=2)
    #axes.errorbar([],[],markerfacecolor='None',ls='',label=f"{leg1}\n{leg2}")
    #axes.set_xscale('log')
    #axes.set_yscale('log')    
    #axes.legend(loc='best', frameon=True)
    ## setting common x and y axes lables.
    #fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
    #fig.text(0.04, 0.5, '$\Delta\Sigma$ [hM$_\odot$/pc$^2]$', va='center', rotation='vertical')#, fontsize=16)
    ## for common title to all the subplots
    ##plt.suptitle(r'samples used for weaklensing singal calculation (NYU_VAGC,Zehavi,Niladri)')
    ##plt.suptitle('weak lensing signal around NYU-VAGC galaxy sample dr72safe%d in HSC field of sources.\n\n%s, %s'%(tag,leg1,leg2))#, fontsize=15)
    ##plt.show()
    #plt.savefig(f"{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
   
 

