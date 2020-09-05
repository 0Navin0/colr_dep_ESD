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

def plot_esd(rbin, method, single):

    #return a list of pipeline outputs with same rbin
    def pipelineoutput_rbin(rbin):
        return glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/*/run*rbin{rbin}*")

    def pipelineoutput_ifrandom_rbin_pofz(rbin,pofz,ifrandom="noRandom"):
        #in future if you have more than one files for one rbin,pofz...may be say for - lens galaxies from another survey 'type'...
        #...(here: nyu-vagc) then throw all those cases along with 'run*' in the below code. 
        return glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/*/run*{ifrandom}*rbin{rbin}_{pofz}")[0]
    
    #read the theoretical signal
    #fixed dir tree. try not to alter it.---26Aug2020, 03Sept2020
    #base0 = glob(f"/home/navin/git/colr_dep_ESD/{pofz}/run{run}*/signal_dr72safe_bin*")[0]
    #base0 = glob(f"/home/navin/git/colr_dep_ESD/{pofz}/run*{run}*{pofz}")[0] #03Sept2020
    base0 = glob(f"/home/navin/git/colr_dep_ESD/Theoretical_ESD_signal")[0] #03Sept2020
    with open(f'{base0}/esd_{method}_magbinned_colrdep_rbin{rbin}.json','r') as f: 
        dic = json.load(f)
    magbin = np.array(list(dic.keys()))
    binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
    rp = dic['rp']
   
    #Info of lens catalog supplied to the pipeline-->Specific to downloaded nyu_data..No need to change.
    sample_base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/col_dep_filter_fits/post_catalog.dr72safe'
    identifier = [7,8,9,10]
    colour = ['red','blue']
#------------------------------------------------------------------    
    #fixed dir tree. try not to alter it.---26Aug2020
    #open config file
    #no need to store config file for every run. I'll try to extract the required info from the tree structure itself..like below
    #with open(glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*/*config*")[0], "r") as yamlfile:
    #    config = yaml.load(yamlfile)
    #outputdir = glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*")[0] #03Sept2020
         
    #base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/signal_dr72safe'
    #base = glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{pofz}/run{run}*/signal*")[0]+'/signal_dr72safe' #26Aug2020
#-------------now I don't do plotting for specific 'run' but now for all runs in an 'rbin'...so looping over 'outputdir' now.------    

    d = {True:'single', False:'overplot'}

    for outputdir in pipelineoutput_rbin(rbin): 
        #keep all the plots for a given rbin in one place-->easy to compare
        plotdir = f"/home/navin/git/colr_dep_ESD/Obesrved_signal/rbin{rbin}/" + outputdir.split("/")[-1]
        print("\n\n#***********************************************")
        print(f"pipeline outputdir:{outputdir}")
        print(f"ceating {plotdir}")
        call(f"mkdir -p {plotdir}",shell=True)
        ifrandom, pofz = outputdir.split('_')[-3], outputdir.split('_')[-1]
        #observed signal from here...pipeline output 
        base = outputdir + '/signal_dr72safe' #03Sept2020
        #debug step: for which pipeline ouput file?
        print(f"\nifrandom={ifrandom}, pofz={pofz}, rbin={rbin}, method={method}, single={d[single]}")

        for tag in identifier:
            #fig,axes = plt.subplots(nrows=2,ncols=2)
            #for ii,colr in zip(range(len(axes)),colour): #accessing one row at a time. 
            if single==False:
                #for both colr on the same plot 
                fig,axes = plt.subplots(1,1)

            for colr in colour:
                # get info of samples to be plotted from the red/blue filtered fits files
                hdul = fits.open('%s%d.%s.fits'%(sample_base,tag,colr))
                hdr = hdul[0].header
                leg1 = "abs_mag_bin:(%s,%s)"%(hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3])
                leg2 = "z-bin: %s" %hdr['z_bin'].split('-')[0][:4]+'-'+(hdr['z_bin'].split('-')[1][:4])
                leg3 = '%s galaxy sample'%colr
                print(f"#--------------color changed to {colr}-----------")
                print("Infromation of observed signal from weaklens pipeline(Credit:Surhud M.):")

                # get data from weaklensing pipeline output to plot
                if ifrandom == 'random':
                    #random rotated pipeline output
                    errDeltaSigma_rand = np.loadtxt('%s%d_%s.dat'%(base,tag,colr),usecols=(11),unpack=True)
                    print("from 'outputdir' as printed above: 'random'-->%s%d_%s.dat"%(base,tag,colr))
                    # not checking  ~np.isnan(deltasigma_rand) & ~np.isnan(r_rand) & 
                    notnan_r = ~np.isnan(errDeltaSigma_rand) 
                    #not checking for ~(deltasigma_rand<0) because due to random rotations, the signal is expected to become very small. and may even go negative in many cases.
                    #all I care about is the errorbar after doing random rotations.
                    positives_r = ~(errDeltaSigma_rand<0) 

                    #all conditions same like pofz,corrections etc except that 'noRandom' this time.
                    rebase = pipelineoutput_ifrandom_rbin_pofz(rbin,pofz,ifrandom="noRandom") + '/signal_dr72safe' #04Sept2020  
                    deltasigma,r,errDeltaSigma = np.loadtxt('%s%d_%s.dat'%(rebase,tag,colr),usecols=(5,7,11),unpack=True)
                    print("noRandom-->%s%d_%s.dat"%(rebase,tag,colr))
                    print(f"tag={tag}, color={colr}, {leg1}\nNo of nan values in \"{colr}\" deltasigma,r,errDeltaSigma,errDeltaSigma_rand: {np.sum(np.isnan(deltasigma))},{np.sum(np.isnan(r))},{np.sum(np.isnan(errDeltaSigma))}, {np.sum(np.isnan(errDeltaSigma_rand))}")
                    print(f"No of negative points in\"{colr}\" deltasigma,errDeltaSigma,errDeltaSigma_rand: {np.sum(deltasigma<0)},{np.sum(errDeltaSigma<0)},{np.sum(errDeltaSigma_rand<0)}\n" )
                    notnan = ~np.isnan(deltasigma) & ~np.isnan(r) & ~np.isnan(errDeltaSigma) 
                    #let's incldue negative singals too since, the error bar will account for how positive this signal can be!!
                    positives = ~(errDeltaSigma<0) #~(deltasigma<0) 

                    #debug step: 
                    print("Now specific details:")
                    print("\nif the sizes of notnan_r!=notnan and positives_r!=positives; the number of data points that we get to plot further reduce due to random rotations.check for it and discuss!")
                    print("why should random rotations cause negative values of ESD?\n")
                    #print('notnan_r',notnan_r)
                    #print('notnan',notnan)
                    #print('positives_r',positives_r)
                    #print('positives',positives)
                    print("before filter errDeltaSigma_rand, notnan_r, positives_r-->", errDeltaSigma_rand, notnan_r, positives_r)
                    errDeltaSigma_rand = errDeltaSigma_rand[notnan_r & positives_r & notnan & positives]
                    print("after filter errDeltaSigma_rand-->",errDeltaSigma_rand)
 
                    print("\nbefore filter: deltasigma,errDeltaSigma,notnan,positives-->", deltasigma, errDeltaSigma, notnan, positives)
                    deltasigma = deltasigma[notnan_r & positives_r & notnan & positives]
                    errDeltaSigma = errDeltaSigma[notnan_r & positives_r & notnan & positives]                   
                    r = r[notnan_r & positives_r & notnan & positives]
                    print("after filter: deltasigma,errDeltaSigma-->", deltasigma, errDeltaSigma)

                    #debug step: is errDeltaSigma_rand > errDeltaSigma ?
                    #errDeltaSigma not to be used for plotting...use errDeltaSigma_rand
                    print(f"\nafter filter-->is errDeltaSigma_rand > errDeltaSigma ?  {errDeltaSigma_rand > errDeltaSigma}")

                    errDeltaSigma = errDeltaSigma_rand 
                    print("\nshould match the above printed errDeltaSigma_rand", errDeltaSigma)
                else:    
                    deltasigma,r,errDeltaSigma = np.loadtxt('%s%d_%s.dat'%(base,tag,colr),usecols=(5,7,11),unpack=True)
                    print("noRandom-->%s%d_%s.dat"%(base,tag,colr))
                    print(f"\ntag={tag}, color={colr}, {leg1}\nNo of nan values in \"{colr}\" deltasigma,r,errDeltaSigma: {np.sum(np.isnan(deltasigma))},{np.sum(np.isnan(r))},{np.sum(np.isnan(errDeltaSigma))}")
                    print(f"No of negative points in\"{colr}\" deltasigma,errDeltaSigma: {np.sum(deltasigma<0)},{np.sum(errDeltaSigma<0)}\n" )
                    notnan = ~np.isnan(deltasigma) & ~np.isnan(r) & ~np.isnan(errDeltaSigma)
                    #singnal can be negative with some finite errobar which can make the singal meaningful
                    positives = ~(errDeltaSigma<0) #~(deltasigma<0) & 
                    #debug step: 
                    print("Now specific details:")
                    #print('notnan',notnan)
                    #print('positives',positives)
                    
                    print("before filter: deltasigma,errDeltaSigma,notnan,positives-->", deltasigma, errDeltaSigma, notnan, positives)
                    deltasigma=deltasigma[notnan & positives]
                    r = r[notnan & positives]
                    errDeltaSigma = errDeltaSigma[notnan & positives]
                    print("after filter: deltasigma,errDeltaSigma-->", deltasigma, errDeltaSigma)
 
                #get calculated(theoretical) ESD from AUM
                if colr=='red':
                    esd = dic[magbin[1:][binmax==float(hdr['absmmax'])][0]][0]['red']
                else: 
                    esd = dic[magbin[1:][binmax==float(hdr['absmmax'])][0]][1]['blue']

                if single==True:
                    # plotting one colr on one plot
                    fig,axes = plt.subplots(1,1)
                    axes.plot(rp, esd, label=f"Prediction:(by HOD model)Niladri et al.",c=colr) #marker='d', markerfacecolor='white',
                    axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC_data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
                    axes.errorbar([],[],markerfacecolor='None',ls='',label=f"{leg1}\n{leg2}")
                    #axes.scatter([],[],facecolors='None',label="=======")
                    axes.set_ylim(0.01,3000)
                    axes.set_xscale('log')
                    axes.set_yscale('log')    
                    axes.legend(loc='best', frameon=True)
                    plt.title(f"{ifrandom} rotation and {pofz} for the HSC sources in rbin={rbin}")
                    #plt.title(f"{method} params used to get HODs.(Niladri et al.)")
                    # setting common x and y axes lables.
                    fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
                    fig.text(0.04, 0.5, '$\Delta\Sigma$ [hM$_\odot$/pc$^2]$', va='center', rotation='vertical')#, fontsize=16)
                    # for common title to all the subplots
                    plt.savefig(plotdir + f"/{method}_{colr}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
                    print(f"\nsaved plot as: " + plotdir + f"/{method}_{colr}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
                else:
                    # plotting red and blue on the same plot 
                    if colr=='blue':
                        axes.plot(rp, esd, label=f"Prediction:(HOD model)Niladri et al.",c=colr) #marker='d', markerfacecolor='white'
                        axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC_data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
                    else:
                        axes.plot(rp, esd, c=colr) #, marker='d',markerfacecolor='white'
                        axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', fmt='o', linewidth=1.5, capsize=5, capthick=2)
            if single==False:
                plt.title(f"{ifrandom} rotation and {pofz} for the HSC sources in rbin={rbin}")
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
                plt.savefig(plotdir + f"/{method}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
                print(f"\nsaved plot as: " + plotdir + f"/{method}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")

        cmd = "mkdir -p " + plotdir + f"/{d[single]}.{method}"
        call(cmd,shell=True)
        print(f"created "+ plotdir + f"/{d[single]}.{method}" )
        call(f"mv {plotdir}/*png {plotdir}/{d[single]}.{method}",shell=True)
        print(f"moved plots to "+ plotdir + f"/{d[single]}.{method}" )
        call(f"tar -czvf {plotdir}/{d[single]}.{method}.tar.gz {plotdir}/{d[single]}.{method}",shell=True)
        
   
if __name__=="__main__":
    method = ['bestfit','fittingFunc']
    #rbin = 10 #5 #15  #change it as per the config file for weaklens_pipeline---older version.

    for rbin in [3,5,10]:
        for ii in method:
            for jj,tt in zip((True,False),("single","overplot")): 
                #tt unused.
                plot_esd(rbin, method=ii, single=jj) 



"""

At some stage I check in the files outputted from weaklens pipeline:
1)There are some negative values in the esd/deltasigma in the observed signal(which was calulated from weaklens pipeline.)
Do you understan why would esd be negative theoretically? negative values also have some error bars...so they are also plotted.
2) Also there are nan values in deltasigma and errDeltaSigma.
-------> nan values are filtered out.
------->-ve values doesn't cause any problem in plotting since we are not taking log()of these numbers...only scale is in log.
------->But in ratio plots these -ve points are excluded.
"""
