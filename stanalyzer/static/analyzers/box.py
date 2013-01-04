#/usr/bin/env python

# for MDAnalysis
import sys
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import math
import numpy as np

# for server side jobs
import os
import string
import random

# import others
import pprint
import pickle
import re
from datetime import datetime

# for using sqlite3
import sqlite3
import stanalyzer

# for local functions and classes
from stanalyzer import *

#**********************************************************************
# *  This function analyzer the system box through all trajectories
#
#**********************************************************************

#///////////////////////////////////////////////////////////////////////////
# Get job information
# -- Use following codes to make your own function
#///////////////////////////////////////////////////////////////////////////
frame_Unit = 0.005;

exe_file = sys.argv[0];
in_file = sys.argv[1];
#out_file = sys.argv[2];

print 'execution file: {}'.format(exe_file);
print 'input file: {}'.format(in_file);
#print 'output file: {}'.format(out_file);

#print 'Reading pickle...'
fid_in = open(in_file, 'rb');
dic = pickle.load(fid_in);
fid_in.close();

# variables in binary file
para_pkeys      = dic["para_pkeys"];        # gui_parameter primary keys obtained from job submitting
job_pkey        = dic["job_pkey"];          # gui_job primary keys obtained from job submitting
job_title       = dic["job_title"];         # gui_job title
prj_pkey        = dic["pkey"];              # gui_project primary key beloning to current job
prj_title       = dic["ptitle"];            # gui_project title beloning to current job
SESSION_HOME    = dic["session_home"];      # default session home directory: ~/dJango_home/media/user_id
OUTPUT_HOME     = dic["output_home"];       # output directory given by user
PBS_HOME        = dic["pbs_home"];          # directory where PBS script is written (i.e. for using cluster)
ANALYZER_HOME   = dic["analyzer_home"];     # directory where all analyzers are located (i.e. the location of this file)
MEDIA_HOME      = dic["media_home"];        # media directory ~/dJango_home/media
DB_FILE         = dic["dbName"];            # database file name and full location 
base_path       = dic["base_path"];         # base location of input files
path_output     = dic["path_output"];       # the location of output directory
path_python     = dic["path_python"];       # python path to run analyzers
structure_file  = dic["structure_file"];    # full path of structure file (i.e. PDB, PSF, etc)
pbs             = dic["pbs"];               # PBS script for using cluster machine
trajectoryFile = [];                        # the list of trajectory files
for i in range(len(dic["trajectory"])):
    trajectoryFile.append(dic["trajectory"][i]);


# Identifying my parameters
myPara = get_myfunction(exe_file, dic);
fName = myPara[0];      # name of function except '.py'
pInfo = myPara[1];      # parameter information pInfo[0] = "number of parameters"
paras = myPara[2];      # actual parameters paras[0] contains 'the number of parameters'
para_pkey = myPara[3];  # primary key of parameter table contating this analzyer function. 

# display all local variables
#pprint.pprint(dic)

# Updating DB: running Job
# 0 - submit job
# 1 - Running job
# 2 - Error occurred
# 3 - Completed
print "para_pkey: "
print para_pkey
stime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
conn = sqlite3.connect(DB_FILE);
c    = conn.cursor();
for i in range(len(para_pkey)):
    query = """UPDATE gui_parameter SET status = "RUNNING" WHERE id = {}""".format(para_pkey[i]);
    c.execute(query);
    conn.commit();
conn.close();


#///////////////////////////////////////////////////////
# print gui_job and gui_parameter table
#///////////////////////////////////////////////////////
stime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
conn = sqlite3.connect(DB_FILE);
c    = conn.cursor();

print "========= gui_job ==========="
query = "SELECT id, name, proj_id, anaz, status, output, stime, etime FROM gui_job";
c.execute(query);
job = c.fetchall();
print "ID\tTITLE\tPROJ_ID\tANALYZER\tSTATUS\tOUTPUT\tSTART\tEND";
for item in job:
    print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(item[0], item[1], item[2], item[3], item[4], item[5], item[6], item[7]);

print "========= gui_parameter ==========="
query = "SELECT id, job_id, anaz, para, val, status FROM gui_parameter WHERE job_id = {}".format(job_pkey[0]);
print query
c.execute(query);
PR = c.fetchall();
print "ID\tJOB_ID\tANALYZER\tPARAMETER\tVALUE\tSTATUS";
for item in PR:
    print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(item[0], item[1], item[2], item[3], item[4], item[5]);
conn.close();


#///////////////////////////////////////////////////////////////////////////
# Running actual job
#///////////////////////////////////////////////////////////////////////////
out_dir = "{}/{}".format(OUTPUT_HOME, fName[1]);
try:
    # Create output directory
    if not (os.path.isdir(out_dir)):
        print "Creating directory into {}".format(out_dir)
        os.mkdir(out_dir);
    
    
    # -------- Writing input file for web-link
    print "list of PARAMETERS: "
    inFile = '{0}/input.dat'.format(out_dir);
    fid_in = open(inFile, 'w');
    strPara = "Name of Function: {}\n".format(exe_file);
    strPara = strPara + "Total {} parameters \n".format(int(paras[0])+3);
    strPara = strPara + "\t- Base path: {}\n".format(base_path);
    strPara = strPara + "\t- Structure file: {}\n".format(structure_file);
    strPara = strPara + "\t- Trajectory files: \n"
    tmp = "";
    for trj in trajectoryFile:
        tmp = tmp + "\t\t{}\n".format(trj);
        
    strPara = strPara + tmp;
    
    strPara = strPara + "\t- Job specific parameters: \n"
    tmp = "";
    for i in range(len(pInfo)):
        tmp = tmp + "\t\t{}:{}\n".format(pInfo[i], paras[i]);
    strPara = strPara + tmp;
    strPara = strPara + "\nPBS: \n{}\n".format(pbs);
    fid_in.write(strPara);
    fid_in.close();
    
    
    outFile = '{0}/output.dat'.format(out_dir);
    fid_out = open(outFile, 'w')
    fid_out.write("Frame\tX-axis\tY-axis\tZ-axis\tAlpha\tBeta\tGamma\tVolumn\n")
    print '--- Calculating Unicell dimension and volum'
    psf = '{0}{1}'.format(base_path, structure_file);
    print psf;
    cnt = 1;
    for idx in range(len(trajectoryFile)):
        # reading trajectory
        dcd = '{0}{1}'.format(base_path, trajectoryFile[idx]);
        print 'Reading PSF: ' + psf
        print 'Reading DCD: ' + dcd
        u = Universe(psf, dcd);
        print '{0} is done!'.format(idx);
        # read based on frame
        for ts in u.trajectory:
            #t = cnt * frame_Unit;
            t = cnt;
            cnt = cnt + 1;
            # get the class MDAnalysis.coordinates.base.Timestep
            # from MDAnalysis.coordinates.DCD.DCDReader
            ucell = u.trajectory[ts.frame]
            outStr = '{0}'.format(t)
            for j in range(len(ucell.dimensions)):
               outStr = outStr + '\t{0}'.format(ucell.dimensions[j])
            outStr = outStr + '\t{0}\n'.format(ucell.volume)
            fid_out.write(outStr)
            print outStr
    fid_out.close()

    # update gui_parameter & gui_job table when job completed
    etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
    conn = sqlite3.connect(DB_FILE);
    c    = conn.cursor();
    for i in range(len(para_pkey)):
        query = """UPDATE gui_parameter SET status = "COMPLETE" WHERE id = {}""".format(para_pkey[i]);
        print query
        c.execute(query);
        conn.commit();
    
    # update gui_job if every status in gui_parameter are COMPLETE
    query = """SELECT DISTINCT(status) FROM gui_parameter WHERE job_id = {}""".format(job_pkey[0]);
    c.execute(query);
    ST = c.fetchall();
    print "number status = {}".format(len(ST));
    
    if (len(ST) == 1) and (ST[0][0] == "COMPLETE"):
        etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
        query = """UPDATE gui_job SET status = "COMPLETE", etime = "{0}" WHERE id = {1}""".format(etime, job_pkey[0]);
        c.execute(query);
        conn.commit();
        
    conn.close();


#///////////////////////////////////////////////////////////////////////////
# Finalizing  job
# -- Use following codes to make your own function
#///////////////////////////////////////////////////////////////////////////
except:
    # update gui_parameter & gui_job table when job failed 
    etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
    conn = sqlite3.connect(DB_FILE);
    c    = conn.cursor();
    for i in range(len(para_pkey)):
        query = """UPDATE gui_parameter SET status = "FAILED" WHERE id = {}""".format(para_pkey[i]);
        c.execute(query);
        conn.commit();
    query = """UPDATE gui_job SET status = "INTERRUPTED" WHERE id = {1}""".format(job_pkey[0]);
    c.execute(query);
    conn.commit();
    conn.close();




#///////////////////////////////////////////////////////
# print gui_job and gui_parameter table
#///////////////////////////////////////////////////////
stime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
conn = sqlite3.connect(DB_FILE);
c    = conn.cursor();

print "========= gui_job ==========="
query = "SELECT id, name, proj_id, anaz, status, output, stime, etime FROM gui_job";
print "ID\tTITLE\tPROJ_ID\tANALYZER\tSTATUS\tOUTPUT\tSTART\tEND";
c.execute(query);
job = c.fetchall();
for item in job:
    print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(item[0], item[1], item[2], item[3], item[4], item[5], item[6], item[7]);

print "========= gui_parameter ==========="
query = "SELECT id, job_id, anaz, para, val, status FROM gui_parameter";
print "ID\tJOB_ID\tANALYZER\tPARAMETER\tVALUE\tSTATUS";
c.execute(query);
PR = c.fetchall();
for item in PR:
    print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(item[0], item[1], item[2], item[3], item[4], item[5]);

conn.close();
