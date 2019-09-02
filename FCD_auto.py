#!/groups/kemi/freja/anaconda3/bin/python
#Made by Andreas Vishart
#Needs the system seperated into fragments and the keywords: iop(3/33=1) pop=full
#--------------------------------------Imports--------------------------------------------
import numpy as np
import copy as copy
import sys
import argparse

#--------------------------------------Parameters-----------------------------------------
#Investigate MOs
Inv_MO = "n"                                                #Investigate MO "y"/"n"                                           

#Set fragment label
FD=1
FA=2
FC=0                                                        #Chromophore (set to zero if it not exists)
FB=0                                                        #Bridge (set to zero if it not exists)

#Number of investigated MOs
num_inv=8                                                 #2*number+1 of investigated MOs 
#Contribution_criteria
contri_crit=0.40
#Weight of frontier orbitals
wFMO=0.05
#Number of best weight contributions of D, A and C system (auto +1 for D and A) 
len_weight=2
#Energy difference criteria
energy_crit=-0.001
#Charge replacement criteria
Charge_crit=0.1

#print the text ("y"/"n" or "yes"/"no")
print_crit="y"

#--------------------------------------Parser---------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('Input', help="Input file is needed.", type=str)
parser.add_argument("-i", "--Inv_MOs", help="See the different localizations of the MOs (yes/no).", type=str)
parser.add_argument("-fa", "--Frag_A", help="The label of the fragment for the acceptor.", type=int)
parser.add_argument("-fd", "--Frag_D", help="The label of the fragment for the donor.", type=int)
parser.add_argument("-fc", "--Frag_C", help="The label of the fragment for the chromophore.", type=int)
parser.add_argument("-fb", "--Frag_B", help="The label of the fragment for the bridge.", type=int)
args = parser.parse_args()
if args.Inv_MOs:
    Inv_MO = args.Inv_MOs
if args.Frag_A:
    FA = args.Frag_A
if args.Frag_D:
    FD = args.Frag_D
if args.Frag_C:
    FC = args.Frag_C
if args.Frag_B:
    FB = args.Frag_B

#--------------------------------------Load File------------------------------------------
try:
    with open(args.Input) as thefile:
        content = thefile.readlines()                   #Load output file
except:
    sys.exit("Error: Missing output file")

#Search file
number = 0
for linje in content:
    if "alpha electrons" in linje:
        AO_number = number-1
        content_AO=content[number-1:number+1]
    if "Molecular Orbital Coefficients" in linje: 
        MO_number = number
        content_MO = content[number:]
    if "*** Overlap" in linje:
        Overlap_number = number
    if "Density Matrix" in linje:
        dens_number = number
    if "Orbital energies and kinetic energies" in linje: 
        energi_number = number
        content_energi = content[number:]
    if "Symbolic Z-matrix:" in linje:
        frag_number = number
    number += 1

#--------------------------------------Functions------------------------------------------
#Constants
def constants():                                                                      
    class con:
        h    = 6.62606957e-34                           # Planck's constant  
        m_e  = 9.10938291e-31                           # electron mass
        m_u  = 1.66053892e-27                           # atomic mass unit
        k_B  = 1.3806488e-23                            # Boltzmann constant
        e    = 1.60217657e-19                           # electron charge
        Har_to_eV = 27.21138602                         # Hartree to eV
        AU_to_Debyee = 2.541765                         # Atomic unit to debyee
        eV_to_J = 1.6021766208*1e-19                    # eV to Joule
    return (con)
con = constants() 

#Find the Number of basis functions (AO)
def num_AO():                                           
    global content_AO
    AO_line = content_AO[0]
    AO_line = [line for line in AO_line.split(" ") if line.isdigit()]
    numbers_AO = int(AO_line[0])
    return numbers_AO

#Download the Molecular orbital coefficients as a list
def MO_coeff(AO,orbital):                               
    global content_MO
    global MO_number
    number_MO = MO_number
    orbital = str(orbital)
    number = number_MO
    number_coef = 0
    for linje in content_MO:
        if "  "+orbital in linje: 
            if number>number_MO and number_coef == 0 and number_MO != 0: 
                nylinje = content[number]
                if "  "+orbital in nylinje[20:]:       #So it is not the atomic orbital
                    number_coef = number    
                    break
        number += 1 
    if number_coef == 0:
        print("Error: MOs do not exist")
        return []
    alpha = list(map(int,filter(None,content[number_coef].split(" "))))
    sojle = alpha.index(int(orbital))
    MOcoef = []
    for i in range(0,AO):
        beta = content[number_coef+3+i][23:]
        beta = list(map(float,filter(None,beta.split(" "))))
        MOcoef.append(beta[sojle])
    return MOcoef

#Overlap integral matrix
def Overlap_Matrix(AO):                                     
    global content
    global Overlap_number
    num_over = Overlap_number+2
    Matrix = np.zeros((AO,AO))
    rows = int(AO/5)
    if AO % 5 > 0:
        rows += 1
    k = 0
    for r in range(0,rows):
        for i in range(0,AO-5*r):
            beta = content[num_over+i+k].replace("D","e")
            beta = list(map(float,filter(None,beta.split(" "))))
            beta = beta[1:]
            for j in range(0,len(beta)):
                Matrix[i+5*r,j+5*r] = beta[j]
                Matrix[j+5*r,i+5*r] = beta[j]
        k += AO+1-5*r
    return Matrix

#Charge difference and transition
def Charge_prop(AO,A_MO,D_MO,Frag_A,Frag_D,over_mat):                
    MOA=MO_coeff(AO,str(A_MO)); MOD=MO_coeff(AO,str(D_MO))
    delta_A=0; trans_A=0; delta_D=0; trans_D=0
    for i in Frag_A:
        i = i-1
        for j in range(AO):
            delta_A+=(MOA[i]*MOA[j]-MOD[i]*MOD[j])*over_mat[i,j]
            trans_A+=(MOA[i]*MOD[j]+MOD[i]*MOA[j])*over_mat[i,j]
    for i in Frag_D:
        i = i-1
        for j in range(AO):
            delta_D+=(MOA[i]*MOA[j]-MOD[i]*MOD[j])*over_mat[i,j]
            trans_D+=(MOA[i]*MOD[j]+MOD[i]*MOA[j])*over_mat[i,j]
    DeltaQ=delta_A-delta_D
    transQ=0.5*(trans_A-trans_D)
    return DeltaQ,transQ

#Energy of a single MO
def Energy(AO,MO):                              
    global content_energi
    global energi_number
    MO = str(MO)
    number_energi = energi_number
    number = number_energi
    number_coef=0
    for linje in content_energi:
        if " "+MO+"   " in linje: 
            if number>number_energi and number_coef == 0 and number_energi != 0: 
                number_coef = number
        if number_coef != 0:
            break    
        number+= 1
    if number_coef == 0:
        print("Error: Energies can not be found")
        return 0
    E_MO = list(filter(None,content[number_coef].split(" ")))
    E_MO = (float(E_MO[2]))*con.Har_to_eV
    return E_MO

#Difference in energy
def Delta_E(AO,A_MO,D_MO):                              
    global content_energi
    global energi_number
    A_MO = str(A_MO)
    D_MO = str(D_MO)
    number_energi = energi_number
    number = number_energi
    number_coef=0; number_coef2=0
    for linje in content_energi:
        if " "+A_MO+"   " in linje: 
            if number>number_energi and number_coef == 0 and number_energi != 0: 
                number_coef = number
        if " "+D_MO+"   " in linje: 
            if number>number_energi and number_coef2 == 0 and number_energi != 0: 
                number_coef2 = number
        if number_coef != 0 and number_coef2 != 0:
            break    
        number+= 1
    if number_coef == 0 or number_coef2 == 0:
        print("Error: Energies can not be found")
        return None
    E_A_MO = list(filter(None,content[number_coef].split(" ")))
    E_D_MO = list(filter(None,content[number_coef2].split(" ")))
    dE = (float(E_A_MO[2])-float(E_D_MO[2]))*con.Har_to_eV
    return dE

#Seperate Fragments by setting the fragment label
def Frag_set(AO,FD,FA,FC=0):               
    global frag_number
    number_frag = frag_number
    fragA=[]; fragD=[]; fragB=[]
    for i in range(1,AO):
        con = content[number_frag+1+i]
        con = list(filter(None,con.split(" ")))
        if "fragment="+str(FA) in con[0].lower():
            fragA.append(i)
        elif "fragment="+str(FD) in con[0].lower():
            fragD.append(i)
        elif "fragment="+str(FC) in con[0].lower():
            fragB.append(i)
        elif len(con)==1:
            break
    if len(fragA) == 0 and len(fragD) == 0 and len(fragB) == 0:
        sys.exit("Error: No Fragments in file")
    Fraglist = Frag_extract(AO,fragA,fragD,fragB)
    return Fraglist[0],Fraglist[1],Fraglist[2]             #A,D,B

#Seperate Fragments by giving the atoms in fragments
def Frag_extract(AO,A_frag,D_frag,B_frag):              
    global content
    global MO_number
    number_MO = MO_number
    Frag_A=[]; Frag_D=[]; Frag_C=[]
    number_coef = 0
    for i in range(0,AO):
        con = content[number_MO+4+i]
        con = con[:21]
        con = list(filter(None,con.split(" ")))
        if len(con) >= 4:
            if int(con[1]) in D_frag:
                number_coef = 1
                Frag_D.append(int(con[0]))
            elif int(con[1]) in A_frag:
                number_coef = 2
                Frag_A.append(int(con[0]))
            elif int(con[1]) in B_frag:
                number_coef = 3
                Frag_C.append(int(con[0]))
            else:
                number_coef = 0
        elif len(con) <= 3:
            if number_coef == 1:
                Frag_D.append(int(con[0]))
            elif number_coef == 2:
                Frag_A.append(int(con[0])) 
            elif number_coef == 3:
                Frag_C.append(int(con[0])) 
            else:
                number_coef = 0
    return Frag_A,Frag_D,Frag_C   

#Investigate different MO in percentage
def perc_MO(AO,number_inter,MO_list,contri_A,contri_D,contri_C):
    for j in range(len(MO_list)):
        if j<number_inter:
            orb_txt="HOMO"
            orb_num=str(j-(number_inter-1))
        else:
            orb_txt="LUMO+"
            orb_num=str(j-(number_inter))
        if len(contri_C)>0:
            print(str(MO_list[j])+" "+orb_txt+orb_num+" :   E="+"{0:.4f}".format(round(Energy(AO,MO_list[j]),4))+ " ;  A: "+str(round(contri_A[j]*100,2))+"%  ;  D: "+str(round(contri_D[j]*100,2))+"%  ;  C: "+str(round(contri_C[j]*100,2))+"%")
        else:
            print(str(MO_list[j])+" "+orb_txt+orb_num+" :   E="+"{0:.4f}".format(round(Energy(AO,MO_list[j]),4))+ " ;  A: "+str(round(contri_A[j]*100,2))+"%  ;  D: "+str(round(contri_D[j]*100,2))+"%")


#MO weight function calculator
def MO_weight_function(MOS_weightest2,wFMO,FMO,MO,Frag,vir_mo_num,contri,contri_crit):
    if contri[-1]>=contri_crit:
        weight2=contri[-1]-wFMO*abs(FMO-MO)
        if weight2>MOS_weightest2[Frag+vir_mo_num+"_weight"][-1]:
            MOS_weightest2[Frag+vir_mo_num+"_weight"][-1]=weight2
            MOS_weightest2[Frag+vir_mo_num][-1]=MO
            MOS_weight_sorted=sorted(zip(MOS_weightest2[Frag+vir_mo_num+"_weight"],MOS_weightest2[Frag+vir_mo_num]),reverse=True)
            MOS_weightest2[Frag+vir_mo_num+"_weight"],MOS_weightest2[Frag+vir_mo_num]=[list(k) for k in zip(*MOS_weight_sorted)]
    return MOS_weightest2


#Calculate the electronic couplings automaticily 
def FCD_auto(Inv_MO,FD,FA,FC,FB,number_inter,contri_crit,wFMO,len_weight,energy_crit,Charge_crit,print_crit):                        
    global content_AO
    AO=num_AO()
    #Dictionary of couplings properties
    ECs={"Coupling_Name":[],"Coupling":[],"A_MO":[],"D_MO":[],"Contribution_D":[],"Contribution_A":[],"Weight_Average":[],"Charge_replacement":[],"Energy_difference":[],"Print_txt":[]}
    #Overlap matrix
    over_mat=Overlap_Matrix(AO)
    #Find the LUMO orbital
    LUMO=int(list(filter(None,content_AO[1].split(" ")))[0])+1
    #Make the list of MO of interest
    MO_list = [(LUMO-number_inter+j) for j in range(0,2*number_inter+1)]
    if FB!=0 and FC==0:
        Frag_A,Frag_D,Frag_C = Frag_set(AO,FD,FA,FB)
    else:
        Frag_A,Frag_D,Frag_C = Frag_set(AO,FD,FA,FC)
    contri_A=[]; contri_D=[]; contri_C=[]
    len_Frag_C=len(Frag_C)
    #How many couplings must be considered
    if len_Frag_C==0:
        len_weight+=1 ; contri_crit=0.0
    MOS_weightest={"A":[0]*len_weight,"A_weight":[0]*len_weight,"A*":[0]*len_weight,"A*_weight":[0]*len_weight,"D":[0]*len_weight,"D_weight":[0]*len_weight,
    "D*":[0]*len_weight,"D*_weight":[0]*len_weight,"C":[0]*len_weight,"C_weight":[0]*len_weight,"C*":[0]*len_weight,"C*_weight":[0]*len_weight}
    vir_mo_num=""; FMO=LUMO-1
    for j in MO_list:
        MO=np.array(MO_coeff(AO,j))**2
        tot_MO=sum(MO)
        #Calculate the contribution of the fragment for each MO
        contri_A.append(sum([(MO[i-1]) for i in Frag_A])/tot_MO)
        contri_D.append(sum([(MO[i-1]) for i in Frag_D])/tot_MO)
        if j>=LUMO:
            vir_mo_num="*"; FMO=LUMO
        #If a chromophore is used
        if len_Frag_C>0:
            contri_C.append(sum([(MO[i-1]) for i in Frag_C])/tot_MO)
            MOS_weightest=MO_weight_function(MOS_weightest,wFMO,FMO,j,"C",vir_mo_num,contri_C,contri_crit)
        #Calculate the weight
        MOS_weightest=MO_weight_function(MOS_weightest,wFMO,FMO,j,"A",vir_mo_num,contri_A,contri_crit)
        MOS_weightest=MO_weight_function(MOS_weightest,wFMO,FMO,j,"D",vir_mo_num,contri_D,contri_crit)
    if Inv_MO.lower()=="y" or Inv_MO.lower()=="yes":
        perc_MO(AO,number_inter,MO_list,contri_A,contri_D,contri_C)
        return " "
    #The electronic coupling for system without chromophore
    ECs=ECs_different(AO,ECs,MOS_weightest,len_Frag_C,FC,FB,Frag_A,Frag_D,Frag_C,number_inter,energy_crit,Charge_crit,contri_A,contri_D,contri_C,LUMO,over_mat,print_crit)
    return ECs

#Which different electronic couplings to calculate
def ECs_different(AO,ECs,MOS_weightest,len_Frag_C,FC,FB,Frag_A,Frag_D,Frag_C,number_inter,energy_crit,Charge_crit,contri_A,contri_D,contri_C,LUMO,over_mat,print_crit):
    if print_crit.lower()=="y" or print_crit.lower()=="yes":
        print("Coupling ; EC value ; D->A ; contri D ; contri A ; avg weight ; charge replaced ; dE")
    #Donor and Acceptor system
    if len_Frag_C==0:
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"A","D",Frag_A,Frag_D,number_inter,energy_crit,Charge_crit,contri_A,contri_D,LUMO,over_mat,print_crit,"H_DA")
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"A*","D*",Frag_A,Frag_D,number_inter,energy_crit,Charge_crit,contri_A,contri_D,LUMO,over_mat,print_crit,"H_D*A*")
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"D","A*",Frag_D,Frag_A,number_inter,energy_crit,Charge_crit,contri_D,contri_A,LUMO,over_mat,print_crit,"H_A*D")
    #Donor, Acceptor and Chromophore system
    elif len_Frag_C>0 and FC!=0:
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"C","D",Frag_C,Frag_D,number_inter,energy_crit,Charge_crit,contri_C,contri_D,LUMO,over_mat,print_crit,"H_DC")
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"A*","C*",Frag_A,Frag_C,number_inter,energy_crit,Charge_crit,contri_A,contri_C,LUMO,over_mat,print_crit,"H_C*A*")
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"D","A*",Frag_D,Frag_A,number_inter,energy_crit,Charge_crit,contri_D,contri_A,LUMO,over_mat,print_crit,"H_A*D")
    #Donor, Acceptor, and Bridge system
    elif len_Frag_C>0 and FB!=0:
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"A*","D*",Frag_A,Frag_D,number_inter,energy_crit,Charge_crit,contri_A,contri_D,LUMO,over_mat,print_crit,"H_D*A*")
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"C*","D*",Frag_C,Frag_D,number_inter,energy_crit,Charge_crit,contri_C,contri_D,LUMO,over_mat,print_crit,"H_D*B*")
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"A*","C*",Frag_A,Frag_C,number_inter,energy_crit,Charge_crit,contri_A,contri_C,LUMO,over_mat,print_crit,"H_B*A*")
        ECs=pairing_fragment(AO,ECs,MOS_weightest,"D","A*",Frag_D,Frag_A,number_inter,energy_crit,Charge_crit,contri_D,contri_A,LUMO,over_mat,print_crit,"H_A*D")
    return ECs

#Pairing the fragments and then calculate the electronic coupling in the order of average weight
def pairing_fragment(AO,ECs,MOS_weightest,Label2,Label1,Frag_2,Frag_1,number_inter,energy_crit,Charge_crit,contri_2,contri_1,LUMO,over_mat,print_crit,EC_txt):
    avg_weight=[]
    MO_coup=[]
    #Get the avgerage weight and the MOs coupled
    for i in range(len(MOS_weightest[Label1])):
        for j in range(len(MOS_weightest[Label2])):
            avg_weight.append(100*(float(MOS_weightest[Label1+"_weight"][i])+float(MOS_weightest[Label2+"_weight"][j]))/2)
            MO_coup.append([MOS_weightest[Label1][i],MOS_weightest[Label2][j]])
    #Get the order
    order_list=sorted(zip(avg_weight,MO_coup),reverse=True)
    #Sort the avgerage weight from highest to smallest and the filenames
    avg_weight,MO_coup=[list(k) for k in zip(*order_list)]
    #Calculate electronic coupling
    for i in range(len(MO_coup)):
        ECs=E_C_FCD(AO,ECs,MO_coup[i][1],MO_coup[i][0],Frag_2,Frag_1,number_inter,energy_crit,Charge_crit,contri_2,contri_1,LUMO,EC_txt,over_mat,print_crit,avg_weight[i])
    if print_crit.lower()=="y" or print_crit.lower()=="yes":
        print("")
    return ECs


#Energy and Charge replacement check, and then electronic coupling
def E_C_FCD(AO,ECs,A_MO,D_MO,Frag_A,Frag_D,number_inter,energy_crit,Charge_crit,contri_A,contri_D,LUMO,EC_text,over_mat,print_crit,avg_weig):
    #The donor and acceptor MO must be different
    if A_MO!=D_MO:
        if A_MO!=0 and D_MO!=0:
            #Energy change
            deE=Delta_E(AO,A_MO,D_MO)
            #The energy criteria has to be fulfilled
            if deE<energy_crit:
                DelQ,tranQ=Charge_prop(AO,A_MO,D_MO,Frag_A,Frag_D,over_mat)
                char_repl=(DelQ**2+4*tranQ**2)**0.5
                #The charge replacement criteria has to be fulfilled
                if abs(char_repl)>Charge_crit:
                    #Electronic coupling in eV
                    H_DA=abs(tranQ)*deE/char_repl
                    #EC in meV
                    H_DA=H_DA*1000
                    #Add to dictionary
                    ECs["Coupling"].append(H_DA)
                    ECs["Charge_replacement"].append(char_repl)
                    ECs["Coupling_Name"].append(EC_text) ; ECs["A_MO"].append(A_MO) ; ECs["D_MO"].append(D_MO) ; ECs["Energy_difference"].append(deE)
                    ECs["Contribution_A"].append(contri_A[A_MO-(LUMO-number_inter)]*100) ; ECs["Contribution_D"].append(contri_D[D_MO-(LUMO-number_inter)]*100)
                    ECs["Weight_Average"].append(avg_weig)
                    #Print text
                    print_txt=str(EC_text)+": "+str(H_DA)+" meV"+" ; "+str(D_MO)+"->"+str(A_MO)+" ; "+str(round(ECs["Contribution_D"][-1],2))+"% ; "+str(round(ECs["Contribution_A"][-1],2))+"% ; "+str(round(ECs["Weight_Average"][-1],2))+"% ; "+str(round(char_repl,3))+" ; "+str(round(deE,3))
                    ECs["Print_txt"].append(print_txt)
                else:
                    print_txt="Not an EC for "+str(D_MO)+"->"+str(A_MO)+", due to the charge replacement, "+str(round(abs(char_repl),3))+" < "+str(Charge_crit)
            else:
                print_txt="Not an EC for "+str(D_MO)+"->"+str(A_MO)+", due to the energy difference, "+str(round(deE,3))+" > "+str(energy_crit)
            if print_crit.lower()=="y" or print_crit.lower()=="yes":
                print(print_txt)
    return ECs
 
#--------------------------------------Calculations---------------------------------------

if __name__ == "__main__":
    ECs = FCD_auto(Inv_MO,FD,FA,FC,FB,num_inv,contri_crit,wFMO,len_weight,energy_crit,Charge_crit,print_crit)


