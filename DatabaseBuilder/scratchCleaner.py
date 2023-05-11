import sys
import re
import os
import time
import logging

inputFolder=sys.argv[1]
keepWatching=(sys.argv[2] == "watch")
outputFolder=inputFolder


logFile=outputFolder+"watch.log"
logging.basicConfig(filename=logFile,filemode='a', level=logging.DEBUG ,format = '%(asctime)s - %(levelname)s: %(message)s',datefmt = '%m/%d/%Y %I:%M:%S %p')


folders=["COMPLETED","PREEMPTED","FAILED","TIMEOUT"]
for outF in folders:
    newFolder=outputFolder+outF
    if not os.path.exists(newFolder):
        os.makedirs(newFolder)

pattern=re.compile("/scratch/mshokrof/[^ /]*/")



while True:
    logging.info("Waking up")
    files=filter(lambda f: f[-4:] == ".err" , os.listdir(inputFolder) )
    for f in files:
        f=inputFolder+f
        server=f.split("/")[-1].split(".")[1]
        scratchFolder=""
        status="FAILED"
        for line in open(f):
            if scratchFolder == "":
                match=re.search(pattern,line)
                if match:
                    scratchFolder=match[0]

        line=line.strip() 
        if line == "1 of 1 steps (100%) done":
            status="COMPLETED"
        elif "PREEMPTION" in line:
            status="PREEMPTED"
        elif line == "Exiting because a job execution failed. Look above for error message":
            status="FAILED"
        elif "DUE TO TIME LIMIT" in line:
            status="TIMEOUT"
        else:
            status="UNKNOWN"



        if status not in ["UNKNOWN"] and scratchFolder != "":
            logging.info(f + " "+status)
            queue="high2"
            if "bm" in server:
                queue="bmh"
            cleanCommand= "srun -p %s -w %s -t 1:00:00 -c 1  --mem=1G  rm -rf %s"%(queue,server, scratchFolder)
            logging.info(cleanCommand)
            os.system(cleanCommand)
        if status in ["COMPLETED","PREEMPTED","FAILED","TIMEOUT"]:
            outFile=".".join(f.split(".")[:-1])+".out"
            mvCommand= "mv %s %s %s%s"%(f,outFile,outputFolder,status)
            logging.info(mvCommand)
            os.system(mvCommand)
            
    if not keepWatching:
        break
    sleepingTime=15
    logging.info("Cleaning done! going to sleep for %d mins"% sleepingTime)
    time.sleep(sleepingTime*60)
