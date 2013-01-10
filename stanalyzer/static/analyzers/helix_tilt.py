#/usr/bin/env python

#----------------------------------------------------------------------------------------
# Author: Jong Cheol Jeong (korjcjeong@yahoo.com, people.eecs.ku.edu/~jjeong)
# 	  Bioinformatics center, The University of Kansas
#----------------------------------------------------------------------------------------

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
import subprocess

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

print "NAME OF FUNCTION: {}".format(fName);
print "PARAMETER INFO: {}".format(pInfo);
print "PARMETERS: {}".format(paras);
print "PARAMETER PRIMARY KEYS: {}".format(para_pkey);

# display all local variables
#pprint.pprint(dic)

# Updating DB: running Job
# 0 - submit job
# 1 - Running job
# 2 - Error occurred
# 3 - Completed
stime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
conn = sqlite3.connect(DB_FILE);
c    = conn.cursor();
for i in range(len(para_pkey)):
    query = """UPDATE gui_parameter SET status = "RUNNING" WHERE id = {}""".format(para_pkey[i]);
    c.execute(query);
    conn.commit();
conn.close();

out_dir = "{}{}".format(OUTPUT_HOME, fName[1]);
# Create output directory
if not (os.path.isdir(out_dir)):
    #print "Creating directory into {}".format(out_dir)
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
######################################## PLEASE DO NOT MODIFY ABOVE THIST LINE!!!! ############################################









#///////////////////////////////////////////////////////////////////////////
# Running actual job
#///////////////////////////////////////////////////////////////////////////
run = 1;
if run:
    outFile = '{0}/output.dat'.format(out_dir);
    fid_out = open(outFile, 'w')
    fid_out.write("# Frame\tTilt angle\n")
    print '--- Calculating Helix tilt'
    psf = '{0}{1}'.format(base_path, structure_file);
    print psf;
    cnt = 1;
    segid      = paras[1];
    seg_name   = paras[2]    
    st_res     = paras[3];
    ed_res     = paras[4];
    
    # default frame unit
    geo = geomatric();

    # select helix
    selQry = 'segid {0} and resid {1}:{2} and name CA'.format(seg_name, st_res, ed_res);
    
    # calculating tilt using PCA    
    for idx in range(len(trajectoryFile)):
	# reading trajectory
	trj = '{0}/{1}'.format(base_path, trajectoryFile[idx]);
	print 'Reading PSF: ' + psf
	print 'Reading DCD: ' + trj
	u = Universe(psf, trj);
	print '{0} is done!'.format(idx);
	# read based on frame
	for ts in u.trajectory:
            tclock = cnt;
            cnt = cnt + 1;
	    # selecting atoms beloning to the helix
	    selAtoms = u.selectAtoms(selQry);
	    print 'pass: selAtoms'
	    if len(selAtoms) > 1:
		# get principal axis
		p1, pe1, pe2 = selAtoms.principalAxes();
		p2 = np.zeros(3);
		p3 = np.array([0, 0, 1]);
		theta = geo.angle(p1,p2,p3);		# calculating tilt
		outStr = '{0}\t{1}\n'.format(tclock,theta);	# print time and degree
	    else:
		outStr = '[{0}/{1}] does not have CA atoms'.format(ts.frame, len(u.trajectory))
            fid_out.write(outStr)
            print outStr;
    fid_out.close()

    # -------- Drawing graphs
    # Writing Gnuplot script
    outScr = '{0}/gplot.p'.format(out_dir);
    outImg  = '{0}.png'.format(exe_file[:len(exe_file)-3]);
    imgPath = "{0}/{1}".format(out_dir, outImg);
    fid_out = open(outScr, 'w');
    gScript = "set terminal png\n";
    gScript = gScript + "set xlabel 'Frame'\n";
    gScript = gScript + "set ylabel 'Tilt Angle'\n";
    gScript = gScript + "set output '{0}'\n".format(imgPath);
    gScript = gScript + """plot "{0}/output.dat" using 1:2 title "Helix Tilt" with linespoints\n""".format(out_dir);
    fid_out.write(gScript);
    fid_out.close()
    
    # Drawing graph with gnuplot
    subprocess.call(["gnuplot", outScr]);
    
    # gzip all reaults
    outZip = "{0}{1}.tar.gz".format(OUTPUT_HOME, fName[1]);
    subprocess.call(["tar", "czf", outZip, out_dir]);

    # Insert values into gui_outputs
    conn = sqlite3.connect(DB_FILE);
    c    = conn.cursor();
    job_title_func = "{0}: {1}".format(job_title, exe_file[:len(exe_file)-3]);
    query = """INSERT INTO gui_outputs (job_id, name, img, txt, gzip) VALUES ({0}, "{1}", "{2}", "{3}", "{4}")""".format(job_pkey[0], job_title_func, imgPath, outFile, outZip);
    c.execute(query);
    conn.commit();
    conn.close();





######################################## PLEASE DO NOT MODIFY BELOW THIST LINE!!!! ############################################
# update gui_parameter & gui_job table when job completed
    etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
    conn = sqlite3.connect(DB_FILE);
    c    = conn.cursor();
    for i in range(len(para_pkey)):
        query = """UPDATE gui_parameter SET status = "COMPLETE" WHERE id = {0}""".format(para_pkey[i]);
        print query
        c.execute(query);
        conn.commit();
    
    # update gui_job if every status in gui_parameter are COMPLETE
    query = """SELECT DISTINCT(status) FROM gui_parameter WHERE job_id = {0}""".format(job_pkey[0]);
    c.execute(query);
    ST = c.fetchall();
    print "number status = {}".format(len(ST));
    
    if (len(ST) == 1) and (ST[0][0] == "COMPLETE"):
        etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
        query = """UPDATE gui_job SET status = "COMPLETE", etime = "{0}" WHERE id = {1}""".format(etime, job_pkey[0]);
        c.execute(query);
        conn.commit();

        # making tar file
        outZip = "{0}job_{1}.tar.gz".format(OUTPUT_HOME, job_pkey[0]);
        subprocess.call(["tar", "czf", outZip, OUTPUT_HOME]);

        # Inserting compressed tar file for all submitted jobs
        final_title = "[** All JOBs **] {0}".format(job_title);
        query = """INSERT INTO gui_outputs (job_id, name, img, txt, gzip) VALUES ({0}, "{1}", "{2}", "{3}", "{4}")""".format(job_pkey[0], final_title, '', '', outZip);
        c.execute(query);
        conn.commit();
	
    conn.close();

try:
    print "okay!";
#///////////////////////////////////////////////////////////////////////////
# Finalizing  job
# -- Use following codes to make your own function
#///////////////////////////////////////////////////////////////////////////
except:
    # update gui_parameter & gui_job table when job failed
    print "ERROR OCCURRED!"
    etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
    conn = sqlite3.connect(DB_FILE);
    c    = conn.cursor();
    for i in range(len(para_pkey)):
        query = """UPDATE gui_parameter SET status = "FAILED" WHERE id = {0}""".format(para_pkey[i]);
        c.execute(query);
        conn.commit();
    query = """UPDATE gui_job SET status = "INTERRUPTED" WHERE id = {0}""".format(job_pkey[0]);
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

