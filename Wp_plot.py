import matplotlib.pyplot as plt 
import numpy as np
import json
from subprocess import call

def plot_wp(method, f_type, cfac):
    #theoretical calcuated values from AUM
    with open(f'./wp/Wp_{method}_magbinned_colrdep_cfac{cfac}.json','r') as f: 
        dic = json.load(f)
    magbin = np.array(list(dic.keys()))
    binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
    rp = dic['rp']

    mags = list(dic.keys())[1:]

    #get measurement data from Niladri's paper: https://arxiv.org/format/1005.2413 
    b_rp, b_23, b_22, b_21, b_20 = np.genfromtxt('./wp.blue.dat',unpack=True)
    b_e23, b_e22, b_e21, b_e20 = np.genfromtxt('./wp.blue.err.dat',unpack=True)
    
    r_rp, r_23, r_22, r_21, r_20 = np.loadtxt('./wp.red.dat',unpack=True)
    r_e23, r_e22, r_e21, r_e20 = np.loadtxt('./wp.red.err.dat',unpack=True)

    fig,ax = plt.subplots(1,1)
    #mkr_style = ['s','p','D','o']
    #mkr_fclr = ['1','1','None','None']
    line_style = ['solid','--','-.',':']
    ax.plot([],[], mfc='r',ls='', label = "AUM(S.M.)")
    for ii,(mag,ls) in enumerate(zip(mags,line_style)):
        if ii==0: #bin-19--20
            ax.plot(rp[4:], np.array(dic[mag][0]['red'][4:])/10, c='red', ls=ls) #shift by 1.0 dex
            ax.plot(rp[4:], np.array(dic[mag][1]['blue'][4:])/10 ,c='blue', label=mag, ls=ls) #shift by 1.0 dex
            ax.errorbar(b_rp, b_20/10, yerr=b_e20/10, color='deepskyblue', marker='s', ls='', elinewidth=1, capsize=2) #, markersize=8 
            ax.errorbar(r_rp, r_20/10, yerr=r_e20/10, color='tomato', marker='s', ls='', elinewidth=1, capsize=2) 
        elif ii==1: #bin-20--21
            #pass
            ax.plot(rp[4:], np.array(dic[mag][0]['red'][4:])/10**0.5 , c='red', ls=ls) #shift by 0.5 dex
            ax.plot(rp[4:], np.array(dic[mag][1]['blue'][4:])/10**0.5 ,c='blue', label=mag, ls=ls) #shift by 0.5 dex
            ax.errorbar(b_rp, b_21/10**0.5, yerr=b_e21/10**0.5, color='deepskyblue', marker='p', ls='', elinewidth=.5, capsize=2) 
            ax.errorbar(r_rp, r_21/10**0.5, yerr=r_e21/10**0.5, color='tomato', marker='p', ls='', elinewidth=.5, capsize=2) 
        elif ii==2: #bin-21--22
            #pass
            ax.plot(rp[4:], np.array(dic[mag][0]['red'][4:]), c='red', ls=ls)# mfc=fclr
            ax.plot(rp[4:], np.array(dic[mag][1]['blue'][4:]) ,c='blue', label=mag, ls=ls) 
            ax.errorbar(b_rp, b_22, yerr=b_e22, color='deepskyblue', marker='D', fillstyle='none', ls='', elinewidth=1, capsize=2) 
            ax.errorbar(r_rp, r_22, yerr=r_e22, color='tomato', marker='D', fillstyle='none', ls='', elinewidth=1, capsize=2) 
        elif ii==3: #bin-22--23
            #pass
            ax.plot(rp[4:], np.array(dic[mag][0]['red'][4:]), c='red', ls=ls)# mfc=fclr
            ax.plot(rp[4:], np.array(dic[mag][1]['blue'][4:]) ,c='blue', label=mag, ls=ls) 
            ebb=ax.errorbar(b_rp, b_23, yerr=b_e23, color='deepskyblue', marker='o', fillstyle='none', ls='', elinewidth=.5, capsize=2)
            ebb[-1][0].set_linestyle(ls) 
            ebr=ax.errorbar(r_rp, r_23, yerr=r_e23, color='tomato', marker='o', fillstyle='none', ls='', elinewidth=.5, capsize=2) 
            ebr[-1][0].set_linestyle(ls) 
    ax.set_xlim(0.1,40)
    ax.set_ylim(0.05,2e5)    
    ax.set_xscale('log')
    ax.set_yscale('log')    
    ax.legend(loc='best', frameon=True)
    fig.text(0.5, 0.02, r"$r_p [h^{-1}$Mpc]", ha='center')#, fontsize=16)
    fig.text(0.04, 0.5, r"$W_p [h^{-1}$Mpc]", va='center', rotation='vertical')#, fontsize=16)
    # for common title to all the subplots
    ax.set_title(f"cfac = {cfac}")
    call(f"mkdir -p ./wp",shell=True)

    ##used to parse json files exported from WebPlotDigitizer
    #with open("wpd_project.json","r") as f:
    #    d = json.load(f)
    #for i in range(len(d['datasetColl'])):
    #    data_size = len(d['datasetColl'][i]['data'])
    #    x = np.zeros(data_size)
    #    y = np.zeros(data_size)
    #    for j in range(data_size-1):
    #        x[j] = d['datasetColl'][i]['data'][j]['value'][0]
    #        y[j] = d['datasetColl'][i]['data'][j]['value'][1]
    plt.savefig(f"./wp/{method}_Wp_allbins_cfac{cfac}.{f_type}")
            
if __name__=="__main__":
    for key in ['bestfit','fittingFunc']:
        plot_wp(key,f_type='png',cfac=1.0)
            
