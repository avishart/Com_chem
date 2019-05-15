#Made by Andreas L. Vishart
#Just use the script in a python enviroment with pandas and subprocess
#Remember to give the path and name of the file saved
#------------Packages--------------
import numpy as np
import os
import sys
#import commands as com
import pandas as pd
import subprocess as sub
from subprocess import Popen,PIPE,check_output

#------------Parameters--------------
#File used with the path in front
thefile="/home/avishart/Overview_of_computational_jobs.csv"
#Hostname
host="avishart"
#Show the jobid
jobid_command="scontrol show jobid "
#File ending
end=[".out",".log"]
#Column names
ColumnTitle=["Jobid","Filename","Description","Method","Status","Path","Done Date"]
#Clear for error and normal terminations in the file (y/n)
clear_file="n"


#------------Program--------------

#Make or load the overview file
try:
    #Dictonary of the file
    content=pd.read_csv(str(thefile))
except:
    #Dictonary of the file
    dic={"Jobid":pd.Series([]),"Filename":pd.Series([]),"Description":pd.Series([]),"Method":pd.Series([]),"Status":pd.Series([]),"Path":pd.Series([]),"Done Date":pd.Series([])}
    content=pd.DataFrame(dic,columns=ColumnTitle)
    content.to_csv(thefile,columns=ColumnTitle)

#Get the queue info for user
output_squeue=str(check_output("squeue -u "+host,shell=True)).split("\\n")
output_squeue_list=[]
for line in output_squeue:
    output_squeue_list.append(list(filter(None,line.split(' '))))
output_squeue_list=output_squeue_list[1:-1]
#Get the jobid number
jobid_num=[str(output_squeue_list[i][0]) for i in range(len(output_squeue_list))]

#Check if some jobs are running or added
for job in jobid_num:
    new_file_test=0
    for job_known in content["Jobid"]:
        if str(job)==str(job_known):
            new_file_test=1
    #Add new jobs to the file
    if new_file_test==0:
        job_detail=str(check_output(jobid_command+job,shell=True)).split("\\n")
        for line in job_detail:
            if "Name=" in line:
                job_name=line.split("=")[-1]
            elif "WorkDir=" in line:
                job_path=line.split("=")[-1]+"/"
        content2=pd.DataFrame({"Jobid":pd.Series([str(job)]),"Filename":pd.Series([str(job_name)]),"Description":pd.Series(["-"]),"Method":pd.Series(["-"]),"Status":pd.Series(["Running"]),"Path":pd.Series([str(job_path)]),"Done Date":pd.Series(["-"])})
        content=content.append(content2,ignore_index=True,sort=False)
        print(str(job_name)+" added to the file")

#Check if some calculations are done 
for index,row in content.iterrows():
    if "Running" in row["Status"]:
        if str(row["Jobid"]) not in jobid_num:
            normal_term=0
            for ending in end:
                try:
                    with open(str(row["Path"])+str(row["Filename"])+str(ending),"r") as check_file:
                        content_check=check_file.readlines()
                    content_check=content_check[-10:]
                except:
                    content_check=[]
                for l in content_check:
                    if "Normal termination" in l:
                        normal_term=1
                        date_list=list(filter(None,l.split("at")[-1].split(" ")))
                    elif "The number of warnings for this run is" in l:
                        normal_term=1
                    elif "PROGRAM ENDED AT" in l:
                        date_list=list(filter(None,l.split("AT")[-1].split(" ")))
            if normal_term==1:
                content.loc[index:index,"Status"]="Done"
                date="".join([str(d)+" " for d in date_list])
                content.loc[index:index,"Done Date"]=date
                print(row["Path"]+row["Filename"]+" is done")
            else:
                content.loc[index:index,"Status"]="Error"
                print("Error in "+row["Path"]+row["Filename"])

#Sort the Dataframe after status and reindex
content=content.sort_values(by=["Status","Jobid"])
content=content.reset_index(drop=True)

#Insert and save the file
content.to_csv(thefile,columns=ColumnTitle)

#Print table
stable=input("Show the file (y/n)?: ")
if stable.lower() == "yes" or stable.lower() == "y":
    print(content.reindex(columns=ColumnTitle))

#Remove all error terminated and normal terminated jobs in the file
if clear_file.lower()=="y":
    with open(thefile) as clear_thefile:
        content_clear=clear_thefile.readlines()
    new_content=[content_clear[0]]
    for line in content_clear[1:]:
        if "Running" in line:
            new_content.append(line)
    newfile=open(thefile,"w")
    for line in new_content:
        newfile.write(line)
    newfile.close()



