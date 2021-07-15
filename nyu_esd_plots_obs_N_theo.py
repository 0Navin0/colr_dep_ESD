# an extension of code in "Tractwise_data/nyu-vagc/iditSmaples/result_weaklen_pipeline/plotSignal1.py"  

# variables to play with: rbin, method, single, ifrandom='noRandom', surveyname='nyu_samples', pofz='fullpofz', run='*', ran_itr='', base0(theoretical esd signal location), base, rebase, plotdir 
# the location of json file for esd signal is fixed now...only need to give it 'run' and 'pofz' to access the right file.

#comment and uncomment relevant lines for using chopped HOD at 2 places by searching for 'standard' #22March21

# name of config file should always end with "config" ------> NOT USING ANYMORE insdead...I access the ifrandom and pofz from the directory structure itself..which is anyways coming from the actual config file

import matplotlib.pyplot as plt 
from glob import glob
import numpy as np
from astropy.io import fits
import json
from subprocess import call
import yaml
import pandas as pd

def select_pipelineoutput_by_args(surveyname, pofz, run, ifrandom, ran_itr, rbin, *subdirs):
    """get a weaklens ouput dir by selecting arguments"""
    #if 'run' for random and noRandom files are the same, then only supply 'run' while calling the func polt_esd().
    #otherwise don't supply 'run' while calling plot_esd(), try to reach the required pipeline output dir without 'run'. 
    #if you really need 'run' to get to the desired location, then supply 'run' argument for 'random' and 'noRandom' cases in pipelineoutput_random_noRandom_pair().

    print(f"{len(subdirs)} number of extra subdirs are given.")
    return  glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{surveyname}/{pofz}/{'/'.join(subdirs)+'/' if len(subdirs) else ''}run{run}_{ifrandom}{ran_itr}_rbin{rbin}_{pofz}")[0]

def pipelineoutput_random_noRandom_pair(surveyname, pofz, run, ifrandom, ran_itr, rbin, *subdirs):
    """ouput: output dir location for type: (random,noRandom) or (noRandom,noRandom)."""
    try:
        default = select_pipelineoutput_by_args(surveyname, pofz, run, ifrandom, ran_itr, rbin, *subdirs) 
        #
        other   = select_pipelineoutput_by_args(surveyname, pofz, run, 'noRandom',  ''  , rbin, *subdirs)
    except:
        sys.exit("No such file found in select_pipelineoutput_by_args().")

    return default,other 

def plot_esd(rbin, method, single, ifrandom='noRandom', surveyname='nyu_samples', pofz='fullpofz', run='*', ran_itr=''):
    """supply  run,random,pfz etc to chose which pipeline output dir gets plotted."""
    #surveyname: to extend this code to include other theoretical HODs(here:Niladri's model for nyu data) and its plot.
    # put '*' for all possibilities of random and pfz

    #read the theoretical signal
    #fixed dir tree. try not to alter it.
    #base0 = glob(f"/home/navin/git/colr_dep_ESD/Theoretical_ESD_signal")[0] #standard place to look for signals
    base0 = glob(f"/home/navin/git/colr_dep_ESD/Theoretical_ESD_signal/Using_chopped_HODs_when_NtotIsLessThan_10**-2_for_nyu_gals")[0] #signals from chopped HODs, 16/03/21
    with open(f'{base0}/esd_{method}_magbinned_colrdep_rbin{rbin}.json','r') as f: 
        dic = json.load(f)
    magbin = np.array(list(dic.keys()))
    binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
    rp = dic['rp']
   
    #Info of lens catalog supplied to the pipeline-->Specific to downloaded nyu_data..No need to change.
    sample_base = '/home/navin/Tractwise_data/nyu-vagc/iditSmaples/col_dep_filter_fits/post_catalog.dr72safe'
    identifier = [7,8,9,10]
    colour = ['red','blue']

    #plot type
    d = {True:'single', False:'overplot'}

    #which pipeline output dir to plot
    ifran_outputdir,other_outputdir = pipelineoutput_random_noRandom_pair(surveyname, pofz, run, ifrandom, ran_itr, rbin)
    #debug step
    if ifrandom == 'noRandom':
        if ifran_outputdir!=other_outputdir:
            sys.exit("ifran_outputdir and other_outputdir should have been the same! check whether you are accessing correct dirs?")

    random, pfz = ifran_outputdir.split('_')[-3], ifran_outputdir.split('_')[-1]
    #keep all the plots for a given rbin in one place-->easy to compare
    #plotdir = f"/home/navin/git/colr_dep_ESD/Obesrved_signal/nyu/rbin{rbin}/" + ifran_outputdir.split("/")[-1] #standard place to store the plots
    plotdir = f"/home/navin/git/colr_dep_ESD/Obesrved_signal/nyu/Using_chopped_HODs_when_NtotIsLessThan_10**-2_for_nyu_gals/rbin{rbin}/" + ifran_outputdir.split("/")[-1] #chopped HODs, 16/03/21
    print("\n\n#***********************************************")
    print(f"pipeline output dirs:\nrandom--->{ifran_outputdir}\nnoRandom--->{other_outputdir}\n")
    print(f"creating {plotdir}")
    call(f"mkdir -p {plotdir}",shell=True)

    #pipeline output: get errorbar from here.
    base = ifran_outputdir + '/signal_dr72safe'
    #debug step: for which pipeline-ouput file?
    print(f"\nmethod={method}, singleplot?={d[single]}")
    print(f"Error bars from: {random}-->%s*"%base)

    for tag in identifier:
        if single==False:
            #for both colr on the same plot 
            fig,axes = plt.subplots(1,1)

        for colr in colour:
            # get info of samples to be plotted from the red/blue-filtered fits files
            # getting the info of color and magnitude cut wise place-holder information from the actual lens data.
            hdul = fits.open('%s%d.%s.fits'%(sample_base,tag,colr))
            hdr = hdul[0].header
            leg1 = "abs_mag_bin:[%s,%s]"%(hdr['absmmax'][-8:-3],hdr['absmmin'][-8:-3])
            leg2 = "z-bin:[%s-%s]" %(hdr['z_bin'].split('-')[0][:4],hdr['z_bin'].split('-')[1][:4])
            leg3 = '%s galaxy sample'%colr
            print(f"#--------------color changed to {colr}-----------")
            print("Infromation of observed signal from weaklens pipeline(Author:Surhud M.):")

            # plot weak lensing signals with intrinsic shears removed: get errDeltaSigma_rand from ifrandom=random files.
            if random == 'random' or random == "random"+f"{ran_itr}":
                #random rotated error bars from covariance matrices
                try:
                    dferr = pd.read_csv(f'{ifran_outputdir}/cov_mat/stddev_nyu')
                except:
                    sys.exit("dferr file doesn't not exist.")
                ##kept just for the idea that errorbars in the random files can be negative.---what does that mean?
                ##not checking for ~(deltasigma_rand<0) because due to random rotations, the signal is expected to become very small. and may even go negative in many cases.

                #kept all conditions same like pofz,corrections etc except that 'noRandom' this time.
                #rebase = pipelineoutput_surveyname_rbin_pofz_ifrandom(surveyname,rbin,pfz,ifrandom="noRandom") + '/signal_dr72safe' #04Sept2020  
                rebase = other_outputdir + '/signal_dr72safe' #21March21  
                deltasigma,r,errDeltaSigma = np.loadtxt('%s%d_%s.dat'%(rebase,tag,colr),usecols=(5,7,11),unpack=True)
                print("signal from: noRandom-->%s%d_%s.dat"%(rebase,tag,colr))
                print(f"tag={tag}, color={colr}, {leg1}\nNo of nan values in \"{colr}\" deltasigma,r,errDeltaSigma: {np.sum(np.isnan(deltasigma))},{np.sum(np.isnan(r))},{np.sum(np.isnan(errDeltaSigma))}")
                print(f"No of negative points in \"{colr}\" deltasigma,errDeltaSigma: {np.sum(deltasigma<0)},{np.sum(errDeltaSigma<0)}\n" )
                notnan = ~np.isnan(deltasigma) & ~np.isnan(r) & ~np.isnan(errDeltaSigma) 
                #let's incldue negative singals too since, the error bar will account for how positive this signal can be!!
                #positives = ~(errDeltaSigma<0) #~(deltasigma<0) 

                #debug step: 
                #print("Now specific details:")
                #print("\nif the sizes of notnan_r!=notnan and positives_r!=positives; the number of data points that we get to plot further reduce due to random rotations.check for it and discuss!")
                #print("why should random rotations cause negative values of ESD?\n")
                #print('notnan_r',notnan_r)
                #print('notnan',notnan)
                #print('positives_r',positives_r)
                #print('positives',positives)
                #print("before filter errDeltaSigma_rand, notnan_r, positives_r-->", errDeltaSigma_rand, notnan_r, positives_r)
                #errDeltaSigma_rand = errDeltaSigma_rand[notnan_r & positives_r & notnan & positives]
                #print("after filter errDeltaSigma_rand-->",errDeltaSigma_rand)
 
                #print("\nbefore filter: deltasigma,errDeltaSigma,notnan,positives-->", deltasigma, errDeltaSigma, notnan, positives)
                deltasigma = deltasigma[notnan] #& positives
                #errDeltaSigma = errDeltaSigma[notnan_r & positives_r & notnan & positives]                   
                r = r[notnan] #& positives
                errDeltaSigma = errDeltaSigma[notnan]                   
                print("after filter: deltasigma, noRandom errDeltaSigma-->", deltasigma, errDeltaSigma)

                #debug step: is errDeltaSigma_rand > errDeltaSigma ?
                #errDeltaSigma not to be used for plotting...use errDeltaSigma_rand
                #print(f"\nafter filter-->is errDeltaSigma_rand > errDeltaSigma ?  {errDeltaSigma_rand > errDeltaSigma}")

                errDeltaSigma = dferr["%d_%s"%(tag,colr)].values 
                print("\nerrDeltaSigma_rand", errDeltaSigma)
            else:    
                deltasigma,r,errDeltaSigma = np.loadtxt('%s%d_%s.dat'%(base,tag,colr),usecols=(5,7,11),unpack=True)
                print("signal from: noRandom-->%s%d_%s.dat"%(base,tag,colr))
                print(f"\ntag={tag}, color={colr}, {leg1}\nNo of nan values in \"{colr}\" deltasigma,r,errDeltaSigma: {np.sum(np.isnan(deltasigma))},{np.sum(np.isnan(r))},{np.sum(np.isnan(errDeltaSigma))}")
                print(f"No of negative points in \"{colr}\" deltasigma,errDeltaSigma: {np.sum(deltasigma<0)},{np.sum(errDeltaSigma<0)}\n" )
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
                print("plotting the following errdeltasigma",errDeltaSigma)
                fig,axes = plt.subplots(1,1)
                axes.plot(rp, esd, label=f"Prediction:HOD model(Paul et al.)",c=colr) #marker='d', markerfacecolor='white',
                axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC-wide data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
                axes.errorbar([],[],markerfacecolor='None',ls='',label=f"{leg1}\n{leg2}")
                #axes.scatter([],[],facecolors='None',label="=======")
                axes.set_ylim(0.01,3000)
                axes.set_xscale('log')
                axes.set_yscale('log')    
                axes.legend(loc='best', frameon=True)
                #plt.title(f"{ifrandom} rotation and {pofz} for the HSC sources in rbin={rbin}")
                #plt.title(f"{method} params used to get HODs.(Niladri et al.)")
                # setting common x and y axes lables.
                fig.text(0.5, 0.02, r"$R [h^{-1} Mpc$]", ha='center')#, fontsize=16)
                fig.text(0.04, 0.5, '$\Delta\Sigma$ [hM$_\odot$/pc$^2]$', va='center', rotation='vertical')#, fontsize=16)
                # for common title to all the subplots
                plt.savefig(plotdir + f"/{method}_{colr}_{rbin}_{pofz}_{random}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
                print(f"\nsaved plot as: " + plotdir + f"/{method}_{colr}_{rbin}_{pofz}_{random}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
            else:
                # plotting red and blue on the same plot 
                if colr=='blue':
                    print("plotting the following errdeltasigma--blue",errDeltaSigma)
                    axes.plot(rp, esd, label=f"Prediction:HOD model(Paul et al.)",c=colr) #marker='d', markerfacecolor='white'
                    axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', label='Observed: HSC-wide data', fmt='o', linewidth=1.5, capsize=5, capthick=2)
                else:
                    print("plotting the following errdeltasigma--red",errDeltaSigma)
                    axes.plot(rp, esd, c=colr) #, marker='d',markerfacecolor='white'
                    axes.errorbar(r, deltasigma, yerr=errDeltaSigma,c=colr, marker='o', fmt='o', linewidth=1.5, capsize=5, capthick=2)
        if single==False:
            #plt.title(f"{ifrandom} rotation and {pofz} for the HSC sources in rbin={rbin}")
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
            plt.savefig(plotdir + f"/{method}_{rbin}_{pofz}_{random}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")
            print(f"\nsaved plot as: " + plotdir + f"/{method}_{rbin}_{pofz}_{random}_{hdr['absmmin'][-8:-3],hdr['absmmax'][-8:-3]}.png")

    cmd = "mkdir -p " + plotdir + f"/{d[single]}.{method}"
    call(cmd,shell=True)
    print(f"created "+ plotdir + f"/{d[single]}.{method}" )
    call(f"mv {plotdir}/*png {plotdir}/{d[single]}.{method}",shell=True)
    print(f"moved plots to "+ plotdir + f"/{d[single]}.{method}" )
    call(f"tar -czvf {plotdir}/{d[single]}.{method}.tar.gz {plotdir}/{d[single]}.{method}",shell=True)
        
   
if __name__=="__main__":
    method = ['bestfit','fittingFunc']
    #rbin = 10 #5 #15  #change it as per the config file for weaklens_pipeline---older version.

    for rbin in [5]: #[3,5,10]
        for ii in method:
            for jj,tt in zip((True,False),("single","overplot")): 
                #tt unused.
                plot_esd(rbin, method=ii, single=jj, ifrandom='noRandom', ran_itr='') #with random, you might need to put ran_itr.



"""

At some stage I check in the files outputted from weaklens pipeline:
1)There are some negative values in the esd/deltasigma in the observed signal(which was calulated from weaklens pipeline.)
Do you understan why would esd be negative theoretically? negative values also have some error bars...so they are also plotted.
2) Also there are nan values in deltasigma and errDeltaSigma.
-------> nan values are filtered out.
------->-ve values doesn't cause any problem in plotting since we are not taking log()of these numbers...only scale is in log.
------->But in ratio plots these -ve points are excluded.
"""

'''A failed attempt to write a general plotting class generalising to different surverynames and x and y axis.

class Plot:
    """Plot based on some selective arguments like run#,pofz,random etc"""
    def __init__(self, run=None, ifrandom=None, rbin=None, pofz=None, sample_type=None, base0=None, method=None):
        """If any of these parameters are not given, will be given wildcard assignment '*'.
           plot will be done for files selected by these parameters.
           In this pipeline, run,sample_type are mandatory inputs. 
           If in future, you need to make them non-mandatory then do: 
           >>> self.run == '*' if run==None
           and to make something else a mandatory input just remove the default wildcard 
           '*' value at appropriate places like in def of 'pipelineoutput_rbin' etc

           Turns out there's no need to change the base code for the above case instead do this:
           >>> g = Plot(run='*', sample_type='hsc',pofz='full*', ifrandom='random', rbin=10)
           >>> print(g.pipeout)
        """
        self.run           = run 
        self.ifrandom      = ifrandom 
        self.rbin          = rbin
        self.pofz          = pofz 
        self.base0         = base0
        self.sample_type   = sample_type 
        self.method        = method

        # Even if i pass some additional kwargs during initiation of the class,
        # and those class variables exist in the namespace of the class, they are 
        # not accessible by doing self.some_arg_definedby_the_below_code
        #for key,value in kwargs.items():
        #    exec(f"self.{key}='{value}'")
        #    print("if the extra passed keyword arguments are supposed to be integers then NOTICE: exec makes all the arguments here as strings.")

        if self.run==None:
            self.run = '*'
        if self.ifrandom==None:
            self.ifrandom = '*'
        if self.rbin==None:
            self.rbin = '*'
        if self.pofz==None:
        self.pofz = '*'

        if sample_type == None:
            sys.exit("'sample_type' is a mandatory input. e.g. sample_type='hsc'")

    #get the pipeline output directories to plot the signal from there.
    #08Feb21
    self.pipeout = glob(f"/home/navin/git/weaklens_pipeline_SM_edited/configs_n_signals/{self.sample_type}/{self.pofz}/run{self.run}_{self.ifrandom}_rbin{self.rbin}_{self.pofz}")
    print(self.pipeout)

    #read the theoretical esd signal
    if (self.base0 == None) & (self.sample_type == 'nyu_vagc'):
        
        self.base0 = glob(f"/home/navin/git/colr_dep_ESD/Theoretical_ESD_signal")[0] 
        print(f"theoretical signal from {self.base0}")
        with open(f'{self.base0}/esd_{self.method}_magbinned_colrdep_rbin{self.rbin}.json','r') as f: 
            self.dic = json.load(f)
        magbin = np.array(list(self.dic.keys()))
        binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
        self.rp = self.dic['rp']
'''
