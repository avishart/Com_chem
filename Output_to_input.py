#Made by Andreas L. Vishart
#Give the script the output file (.out or .log) and the new filename
#------------------------Packages------------------------
import sys
import os


#------------------------Parameters----------------------
#Ending of the new filename 
end="FCD"
#Functional
func="cam-b3lyp"
#Basis set
basis="cc-pvdz"
#Memory for SLURM
mem="2GB"
#Number of processors for SLURM
nprocs="1"
#Calculation type
calc="pop=full iop(3/33=1)"

#------------------------Functions-----------------------
def load_xyz(output_file):
    #Load the output file
    with open(sys.argv[1]) as thefile:
        out_content=thefile.readlines()
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
#Load the output file and get xyz coord
xyz_content=load_xyz(sys.argv[1])

#Make the name of the new file
new_filename=sys.argv[2]


#Make the new input file
newfile=open(new_filename,"w")
newfile.write("%nprocshared="+str(nprocs)+"\n%mem="+str(mem)+"\n")
newfile.write("# "+calc+" "+func+"/"+basis+"\n\nTitel\n\n0 1\n")
newfile.write(xyz_content)
newfile.write("\n")
newfile.close()


