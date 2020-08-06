import matplotlib.pyplot as plt 
import numpy as np
import json

def plot_wp(method):
    with open(f'./Wp_{method}_magbinned_colrdep.json','r') as f: 
        dic = json.load(f)
    magbin = np.array(list(dic.keys()))
    binmax = np.array([float(magbin[ii+1][:5]) for ii in range(len(magbin)-1)])
    rp = dic['rp']

    mags = list(dic.keys())[1:]
    
    fig,ax = plt.subplots(1,1)
    #mkr_style = ['s','p','D','o']
    #mkr_fclr = ['1','1','None','None']
    line_style = ['solid','--','-.',':']
    ax.plot([],[], mfc='None',ls='', label = "AUM(S.M.)")
    for ii,(mag,ls) in enumerate(zip(mags,line_style)):
        if ii==0:
            ax.plot(rp[4:], np.array(dic[mag][0]['red'][4:])/10, c='red', ls=ls) 
            ax.plot(rp[4:], np.array(dic[mag][1]['blue'][4:])/10 ,c='blue', label=mag, ls=ls) 
        elif ii==1:
            ax.plot(rp[4:], np.array(dic[mag][0]['red'][4:])/10**0.5, c='red', ls=ls) 
            ax.plot(rp[4:], np.array(dic[mag][1]['blue'][4:])/10**0.5 ,c='blue', label=mag, ls=ls) 
        else:
            ax.plot(rp[4:], np.array(dic[mag][0]['red'][4:]), c='red', ls=ls)# mfc=fclr, #shift by 0.5 dex
            ax.plot(rp[4:], np.array(dic[mag][1]['blue'][4:]) ,c='blue', label=mag, ls=ls) #shift by 0.5 dex
    ax.set_xlim(0.1,40)
    ax.set_ylim(0.05,2e5)    
    ax.set_xscale('log')
    ax.set_yscale('log')    
    ax.legend(loc='best', frameon=True)
    fig.text(0.5, 0.02, r"$r_p [h^{-1}$Mpc]", ha='center')#, fontsize=16)
    fig.text(0.04, 0.5, r"$W_p [h^{-1}$Mpc]", va='center', rotation='vertical')#, fontsize=16)
    # for common title to all the subplots
    plt.savefig(f"{method}_Wp_allbins.png")
            
if __name__=="__main__":
    for key in ['bestfit','fittingFunc']:
        plot_wp(key)
            
