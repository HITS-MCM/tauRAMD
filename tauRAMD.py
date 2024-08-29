"""
#------------------------------------------------------------------------
#   Computation of the drug-target relative residence times from RAMD simulations
#>  \version{version 1.1 (March. 2020)}
#> <c
#>  Copyright (c) 2020
#>  Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#>  Schloss-Wolfsbrunnenweg 35
#>  69118 Heidelberg, Germany
#>
#>  Please send your contact address to get information on updates and
#>  new features to "mcmsoft@h-its.org". Questions will be
#>  answered as soon as possible.
#>
#>  Authors: Daria Kokh  Daria.Kokh@h-its.org
"""

from matplotlib import *
from matplotlib import gridspec
import  pylab as plt
import seaborn as sns
import numpy as np
from scipy.stats import norm
import sys
import os


soft = "Gr"    #   if Gromacs software was used for RAMD simulations;  otherwise define  soft = 'NAMD' 

def printUsage():
    print ('''\
    NAME

         tauRAMD.py     - computation of residence times using bootstrapping from RAMD output generated in Gromacs

    USAGE

          python tauRAMD.py  input_file[s]

          input files must contain a set of lines extracted from the gromacs (or NAMD) output. Each line contains the number of steps executed before dissociation 

          and has the following format:
          for Gromacs:   “XX/YYYY.out:==== RAMD ==== GROMACS will be stopped after 874650 steps.”   
          for NAMD:      "EXIT: XXXXXXX  > LIGAND EXIT EVENT DETECTED"

    OUTPUT

         residence time with the standard deviation computed for each input_file and an image with histogram representation of the bootstrapping output
    ''')
    
def bootstrapp(t, rounds=50000):
    max_shuffle = rounds
    alpha = 0.8
    sub_set = int(alpha*len(t))
    tau_bootstr = []
    for i in range(1,max_shuffle):
        # generate a sub-set
        np.random.shuffle(t)
        t_b = t[:sub_set]
        # find residence time from a sub-stet
        t_b_sorted_50 =(np.sort(t_b)[int(len(t_b)/2.0-0.5)]+np.sort(t_b)[int(len(t_b)/2)])/2.0
        tau_bootstr.append(t_b_sorted_50)
    return(tau_bootstr)

#---- reading data

d_list = []
if len(sys.argv) < 2:
    printUsage()
    sys.exit(1)
else:
    for i in range(1,len(sys.argv)):
        if os.path.isfile(sys.argv[i]):
            print ("Data found:", sys.argv[i])
            d_list.append(sys.argv[i])
        else:
            print ("Data not found:", sys.argv[i])


times_set = []

for t,d in enumerate(d_list):
    with open(d) as f:
        read_data = f.readlines()
    times = []
    for r in read_data:
        try:
            if soft == "NAMD":
                times.append(int(r[r.find("EXIT:")+6:r.find(">")-2]))   # if NAMD  was used to generate RAMD trajectories
            else:
                times.append(int(r[r.find("after")+6:r.find("steps")-1]))   # if Gromacs was used to generate RAMD trajectories
        except:
            print("Input files dont't have correct format")
            printUsage()
            sys.exit()
#    print(times)
    times = np.asarray(times)/500000.
    times_set.append(times)
    print("************ Dissociation times for "+d+" ****************")
    print(times)


fig  = plt.figure(figsize = (2*len(d_list), 7))
gs = gridspec.GridSpec(nrows=3, ncols=len(d_list), wspace=0.1,hspace=0.6)
mue_set = []   
print("\n ==============   Bootstrapping and computation of tau for each replica==================\n")
for t, times in enumerate(times_set):
#--- do bootstrapping ----
    if len(times) > 8:
        for ti in range(len(times),15): times = np.concatenate((times,[55.0]))
        print(times)
        ax0 = fig.add_subplot(gs[0, t])
        ax0.hist(times,bins=int(len(times)/2),cumulative=True,histtype="step",color='k',lw=1)
        plt.title("raw CDF",fontsize=12)
        ax0.set_xlabel('dissociation time [ns]', fontsize=10)
        ax0.plot([min(times), max(times)],[len(times)/2.0,len(times)/2.0], color='red', alpha = 0.5)
        if (t> 0): ax0.set_yticklabels( [])
        bt2 = bootstrapp(times, rounds=50000)
#--- make a plot------
        bins = 6
        ax1 = fig.add_subplot(gs[1, t])
        ax1.hist(x=bt2,bins=bins, alpha=0.8,density=True,histtype="step")
        mu, std = norm.fit(bt2)
        mue_set.append(np.round(mu,1))
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        #xmax = np.round(max(times))
        x = np.linspace(0.8*xmin, xmax, 100)
        p = norm.pdf(x, mu, std)
        ax1.plot(x, p, 'k', linewidth=2)
        ax1.plot([mu,mu],[0, max(p)], color='red', alpha = 0.5)
        ax1.plot([xmin, xmax],[max(p)/2.0,max(p)/2.0], color='red', alpha = 0.5)
        ax1.plot([0.8*xmin, mu],[max(p),max(p)], color='red', linestyle='dashed',alpha = 0.5)        
        ax1.set_xlabel('res. time [ns]', fontsize=10)
        plt.title("tau distribution",fontsize=12)
        ax1.set_yticks([])

        ax2 = fig.add_subplot(gs[2, t])
        xmin = min(times)
        xmax = np.round(max(times))
        tp = np.linspace(xmin*0.5,xmax*1.5,100)
        poisson = 1-np.exp(-tp/mu) #np.cumsum(1-np.exp(-np.linspace(xmin,xmax,10)/mu))
        points=len(times)
        bins = len(times)
        times = np.asarray(times)
        hist, bin_edges = np.histogram(times,bins=bins)
        hist_center = []
        for i,b in enumerate(bin_edges):
            if i > 0: hist_center.append((bin_edges[i-1]+bin_edges[i])/2.0)
        CD = np.cumsum(hist)/np.max(np.cumsum(hist))
        ax2.scatter(np.log10(np.asarray(hist_center)),CD,marker='o')
        ax2.set_xlabel('log(res. time [ns])', fontsize=10)
        ax2.plot(np.log10(tp),poisson,color = 'k')
        ax2.set_ylim(0,1)
        ax2.set_xlim(-1.5,1.5)
     #   ax2.set_xlim(np.round(np.log10(np.asarray(hist_center))[0],1)-0.1,np.log10(xmax*1.5)) #max(np.log10(np.asarray(hist_center))[-1],np.log10(tp)[-1]))
        ax2.set_yticks(np.linspace(0,1,5))
        if (t> 0): ax2.set_yticklabels( [])
        plt.grid(linestyle = '--',linewidth=0.5)
        ax2.plot([np.log10(mu),np.log10(mu)],[0, 1], color='red', alpha = 0.5)
       # p_mu = 1-np.exp(-1.0)
       # cd_mu = np.argwhere(np.asarray(hist_center)> p_mu)[0][0]
       # if cd_mu > 1:
       #     ks = abs(p_mu - (CD[cd_mu]+CD[cd_mu-1])/2)
       # else:
       #     ks = abs(p_mu - CD[cd_mu])
        KS = np.round(np.max(np.abs(1-np.exp(-(np.asarray(hist_center))/mu) - CD)),2)
        plt.title("KS test:"+str(KS),fontsize=12)
        print(" Relative res. time and SD: ",np.round(mu,2), np.round(std,2),"KS test:",KS)
        print("-----------------------------------------------------------------")


#fig.align_labels()
plt.savefig('res_times.png', bbox_inches='tight',dpi=300)
fig  = plt.figure(figsize = (2*len(d_list), 2))
meanpointprops = dict(linestyle='--', linewidth=1.5, color='firebrick')
medianpointprops = dict(linestyle='-', linewidth=2.0, color='orange')
#plt.yticks(np.linspace(0,ymax,10),np.linspace(0,ymax,10))
plt.boxplot(times_set,showmeans=True, meanline=True,meanprops=meanpointprops,medianprops = medianpointprops, bootstrap=5000) #labels = mue_set)
ymin, ymax = plt.ylim()
plt.ylim=(0,ymax)
plt.grid(linestyle = '--',linewidth=0.5)
plt.yticks(np.linspace(0,int(ymax),min(int(ymax)+1,11)), fontsize=9)
plt.ylabel('residence time [ns]', fontsize=10)
plt.title("Residence times for "+str(t+1)+" replicas, mean: "+str(np.round(np.mean(mue_set),2))+"  std: "+str(np.round(np.std(mue_set),2)),fontsize=10)
plt.savefig('res_times_summary.png', bbox_inches='tight',dpi=300)
print("Residence times for "+str(t+1)+" replicas, mean: "+str(np.round(np.mean(mue_set),2))+"  std: "+str(np.round(np.std(mue_set),2)))


