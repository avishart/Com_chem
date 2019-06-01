#Made by Andreas L. Vishart
#Either give the script the output file(s) (.out or .log) and then the new filename (Auto_resub=False)
#or give the script a range of output file, where it check if it is done (Auto_resub=True)
#A new input file is made, so it can be resubmitted
#------------------------Packages------------------------
import sys
import os

#------------------------Parameters----------------------
#Extra convergence criteria to help convergence (True/False)
Convergence_crit=True
#Automatic resubmit files (True/False)
Auto_resub=True

#------------------------Functions-----------------------
#Check if file is converged 
def check_conv(out_content):
    Check=False
    max_force=["NO"];rms_force=["NO"]
    max_dis=["NO"];rms_dis=["NO"]
    for line in out_content:
        if "Maximum Force" in line:
            max_force.append(list(filter(None,line.replace("\n","").split(" ")))[-1])
        elif "RMS     Force" in line:
            rms_force.append(list(filter(None,line.replace("\n","").split(" ")))[-1])
        elif "Maximum Displacement" in line:
            max_dis.append(list(filter(None,line.replace("\n","").split(" ")))[-1])
        elif "RMS     Displacement" in line:
            rms_dis.append(list(filter(None,line.replace("\n","").split(" ")))[-1])
    if "NO" in max_force[-1] or "NO" in rms_force[-1] or "NO" in max_dis[-1] or "NO" in rms_dis[-1]:
        if "Normal termination" in out_content[-1] or "Error termination" in out_content[-10:-1]:
            Check=True
    return Check

#Get the header for the input file
def load_header(out_content,Convergence_crit):
    Gaus=False
    lines_shift=False
    for line in range(len(out_content)):
        if "Gaussian" in out_content[line]:
            if "***************************" in out_content[line-1]:
                Gaus=True
        elif "Initial Parameters" in out_content[line]:
            Gaus=False
        if Gaus==True:
            if "%nprocshared=" in out_content[line]:
                nprocs=out_content[line].replace("\n","").split("=")[-1]
            if "%mem" in out_content[line]:
                mem=out_content[line].replace("\n","").split("=")[-1]
            if "-------------------" in out_content[line]:
                lines_shift=True
            elif lines_shift==True:
                if "#" in out_content[line]:
                    calc=out_content[line].replace("\n","").split("#")[-1]
            if "Charge = " in out_content[line]:
                newline=list(filter(None,out_content[line].replace("\n","").split(" ")))
                charge=newline[2].replace("\n","")
                multi=newline[-1].replace("\n","")
    if Convergence_crit==True:
        if "int=" not in calc:
            calc=calc+" int=grid=Ultrafine"
        if "scf=" not in calc:
            calc=calc+" scf=(maxcycles=2048)" 
    return nprocs,mem,calc,charge,multi

#Get the xyz data for the last state
def load_xyz(out_content):
    frag=False
    frag_content=[]
    out=False
    for line in out_content:
        #Get fragments for each element
        if "Charge" in line and "Multiplicity" in line:
            frag=True
        elif "Initial Parameters" in line:
            frag=False
        if frag:
            newline=list(filter(None,line.replace("\n","").split(" ")))
            if len(newline)==4:
                frag_content.append(newline[0])
        #Get the final coordinates
        if "Standard orientation:" in line:
            out=True
            l=1
            out_xyz=[]
        elif "Rotational constants" in line:
            out=False
        elif out:
            newline=list(filter(None,line.replace("\n","").split(" ")))
            if newline[0]==str(l):
                out_xyz.append(["{0:.8f}".format(round(float(newline[3]),8)),"{0:.8f}".format(round(float(newline[4]),8)),"{0:.8f}".format(round(float(newline[5]),8))])
                l+=1
    #If elements and coordinates are loaded coorectly then make the xyz input
    xyz_txt=""
    if len(frag_content)==len(out_xyz):
        for i in range(len(frag_content)):
            space=32-len(str(frag_content[i]))-len(str(out_xyz[i][0]))
            xyz_txt+=" "+str(frag_content[i])+" "*space
            for j in range(2):
                space=16-len(str(out_xyz[i][j+1]))
                xyz_txt+=str(out_xyz[i][j])+" "*space
            xyz_txt+=str(out_xyz[i][2])+"\n"
    return xyz_txt




#------------------------Program-------------------------
if Auto_resub==False:
    #Load the output file
    with open(sys.argv[1]) as thefile:
        out_content=thefile.readlines()
    #Get the header for the input file
    nprocs,mem,calc,charge,multi=load_header(out_content,Convergence_crit)
    #Get xyz coord
    xyz_content=load_xyz(out_content)
    #Make the name of the new file
    new_filename=sys.argv[2]
    #Make the new input file
    newfile=open(new_filename,"w")
    newfile.write("%nprocshared="+str(nprocs)+"\n%mem="+str(mem)+"\n")
    newfile.write("#"+calc+"\n\nTitel\n\n"+charge+" "+multi+"\n")
    newfile.write(xyz_content)
    newfile.write("\n")
    newfile.close()

elif Auto_resub==True:
    for outfile in sys.argv[1:]:
        #Load the output file
        with open(outfile) as thefile:
            out_content=thefile.readlines()
        if check_conv(out_content)==True:
            #Get the header for the input file
            nprocs,mem,calc,charge,multi=load_header(out_content,Convergence_crit)
            #Get xyz coord
            xyz_content=load_xyz(out_content)
            #Make the name of the new file
            new_filename=outfile[:-4]+"_new.com"
            #Make the new input file
            newfile=open(new_filename,"w")
            newfile.write("%nprocshared="+str(nprocs)+"\n%mem="+str(mem)+"\n")
            newfile.write("#"+calc+"\n\nTitel\n\n"+charge+" "+multi+"\n")
            newfile.write(xyz_content)
            newfile.write("\n")
            newfile.close()
            print(new_filename)

