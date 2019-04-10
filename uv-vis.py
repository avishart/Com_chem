#Made by Andreas Vishart
#Made by following the exercise 5 from MolStat
# ":" is used to generate subplots, so data after ":" is in same subplot

#------------Package---------------------------
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt

#------------Parameters---------------------------
#Wavelength interval with set values, (from,to,number of points)
wavelength_interval_set=np.linspace(200,700,1002)
#Search word in the file
search="Excited State"
#Full width half maximum [cm^-1]
sigma=3226.00
#Title (if only "" then the name of the first file)
title="test"
#Filetype of plot
end=".pdf"
#Save to file (y/n)
save_file="y"
#Activate oscillator strength (y/n)
osc_act="y"


#------------Constants in SI-units----------------
#Avogadros Number
NA=6.022140857e23
#Elemental charge
e=1.6021766208e-19
#Electron mass
me=9.10938356e-31
#Vacuum permittivity
eps0=8.8541878176e-12
#Speed of light
c=299792458

#Constant k 
k=NA*e**2*np.sqrt(np.log(2)/np.pi)/(2*np.log(10)*me*eps0*c**2)
#print(k)
k2=k/sigma


#------------Functions------------------------------
#Extract excitation data from file
def extract(filename,search):
    try:
        with open(filename) as thefile:
            content=thefile.readlines()
    except:
        print("No file given!")
    #Find the wavelength & oscillator strength of the transitions 
    exi_state=[]
    wave_trans=[]
    osc_trans=[]
    for line in content:
        if search in line:
            state=list(filter(None,line.split(" ")))
            exi_state.append(state)
            wave_trans.append(float(state[6]))
            osc_trans.append(float(state[8][2:]))
    return exi_state,wave_trans,osc_trans


#Calculate the molar absorption coefficients
def epsilon(wavelength_interval,wave_trans,osc_trans,k2,sigma):
    epsi=[]
    for w in wavelength_interval:
        epsi_inter=0
        for j in range(len(wave_trans)):
            epsi_inter+=osc_trans[j]*np.exp(-4*np.log(2)*((1/w-1/wave_trans[j])/(sigma*1e-7))**2)
        epsi.append(k2*epsi_inter)
    return epsi

#Show Uv-Vis plot with oscillator strength
def plot_osc(wavelength_interval,epsi,wave_trans,osc_trans,end):
    fig,ax1=plt.subplots(figsize=(8,6))
    ax2=ax1.twinx()
    ax1.set_xlabel("Wavelength / [nm]")
    ax1.set_ylabel(r"$\epsilon$ / [M$^{-1}$ cm$^{-1}$]",color="b")
    ax1.plot(wavelength_interval,epsi,"b-")
    if osc_act.lower()=="y" or osc_act.lower()=="yes":
        ax2.set_ylabel("Oscillator strength",color="r")
        ax2.bar(wave_trans,osc_trans,width=1,color="r")
    ax1.set_xlim(wavelength_interval[0],wavelength_interval[-1])
    ax1.set_ylim(0,max(epsi))
    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax2.set_ylim(0,1)
    if save_file.lower()=="y" or save_file.lower()=="yes":
        plt.savefig(title,bbox_inches='tight')
    plt.show()
    plt.close()
    return

#Show Uv-Vis plot with oscillator strength for only one data set
def plot_osc_1data(filename,search,end):
    #Extract excitation data from file 1
    exi_state,wave_trans,osc_trans=extract(filename,search)
    #Wavelength interval with automatic values
    wavelength_interval=np.linspace(wave_trans[-1],wave_trans[0]+100,(wave_trans[0]-wave_trans[-1]+1)*3)
    #Calculate the molar absorption coefficitent
    epsi=epsilon(wavelength_interval,wave_trans,osc_trans,k2,sigma)
    #Plot Uv-Vis with oscillator strength
    plot_osc(wavelength_interval,epsi,wave_trans,osc_trans,end)
    return 

#Show Uv-Vis plot of multiply data set
def plot_multi(search,files,wavelength_interval_set,end):
    num_sub=sys.argv.count(":")
    fig,ax=plt.subplots(num_sub+1,sharex=True,sharey=True,figsize=(8,6+2*num_sub))
    wave_list=[]
    epsi_list=[]
    num_sub_now=0
    if num_sub==0:
        ax=[ax]
    ax_y_copy={}
    if osc_act.lower()=="y" or osc_act.lower()=="yes":
        ax_y_copy[num_sub_now]=ax[num_sub_now].twinx()
        ax_y_copy[num_sub_now].set_ylabel("Oscillator strength")
    for i in range(1,len(sys.argv)):
        if sys.argv[i]!=":":
            exi_state1,wave_trans1,osc_trans1=extract(sys.argv[i],search)
            epsi1=epsilon(wavelength_interval_set,wave_trans1,osc_trans1,k2,sigma)
            epsi_list+=epsi1
            ax[num_sub_now].plot(wavelength_interval_set,epsi1,label=sys.argv[i])
            ax[num_sub_now].legend(loc=0)
            ax[num_sub_now].set_ylim(0,max(epsi_list)*1.05)
            ax[num_sub_now].set_xlabel("Wavelength / [nm]")
            ax[num_sub_now].grid(True)
            if osc_act.lower()=="y" or osc_act.lower()=="yes":
                ax_y_copy[num_sub_now].bar(wave_trans1,osc_trans1,width=1)
        elif sys.argv[i]==":":
            num_sub_now+=1
            if osc_act.lower()=="y" or osc_act.lower()=="yes":
                ax_y_copy[num_sub_now]=ax[num_sub_now].twinx()
                ax_y_copy[num_sub_now].set_ylabel("Oscillator strength")
    fig.subplots_adjust(hspace=0.2)
    fig.text(0.04, 0.5, r"$\epsilon$ / [M$^{-1}$ cm$^{-1}$]", va='center', rotation='vertical')
    ax[0].set_xlim(wavelength_interval_set[0],wavelength_interval_set[-1])
    ax[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
    if save_file.lower()=="y" or save_file.lower()=="yes":
        plt.savefig(title,bbox_inches='tight')
    plt.show()
    plt.close()
    return


#------------Program------------------------------


#Only execute if this is the main script
if __name__ == "__main__":
    #Name of plot
    if len(title)==0:
        title=sys.argv[1][:-4]+end
    else:
        title=title+end
    #Plot either one uv-vis with oscillator strength or multi uv-vis without
    if len(sys.argv)<=2:
        #Uv-Vis plot with oscillator strength for only one data set
        plot_osc_1data(sys.argv[1],search,end)
    else:
        #Uv-Vis plot for multi data set
        plot_multi(search,sys.argv,wavelength_interval_set,end)


