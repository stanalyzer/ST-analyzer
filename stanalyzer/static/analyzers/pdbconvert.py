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
from collections import defaultdict	# for initializing dictionary


# for using sqlite3
import sqlite3
import stanalyzer
import lipidArea

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
exe_file = sys.argv[0];
in_file = sys.argv[1];
ST_para_idx = int(sys.argv[2]);         # parameter index for multiple jobs in a same form
fid_in = open(in_file, 'rb');
dic = pickle.load(fid_in);
fid_in.close();

# variables in binary file
ST_para_pkeys      = dic["para_pkeys"];        # gui_parameter primary keys obtained from job submitting
ST_job_pkey        = dic["job_pkey"];          # gui_job primary keys obtained from job submitting
ST_job_title       = dic["job_title"];         # gui_job title
ST_prj_pkey        = dic["pkey"];              # gui_project primary key beloning to current job
ST_prj_title       = dic["ptitle"];            # gui_project title beloning to current job
ST_SESSION_HOME    = dic["session_home"];      # default session home directory: ~/dJango_home/media/user_id
ST_OUTPUT_HOME     = dic["output_home"];       # output directory given by user
ST_PBS_HOME        = dic["pbs_home"];          # directory where PBS script is written (i.e. for using cluster)
ST_ANALYZER_HOME   = dic["analyzer_home"];     # directory where all analyzers are located (i.e. the location of this file)
ST_MEDIA_HOME      = dic["media_home"];        # media directory ~/dJango_home/media
ST_DB_FILE         = dic["dbName"];            # database file name and full location 
ST_base_path       = dic["base_path"];         # base location of input files
ST_path_output     = dic["path_output"];       # the location of output directory
ST_path_python     = dic["path_python"];       # python path to run analyzers
ST_structure_file  = dic["structure_file"];    # full path of structure file (i.e. PDB, PSF, etc)
pdb_file  	   = dic["pdb_file"];          # full path of structure file (i.e. PDB)
ST_pbs             = dic["pbs"];               # PBS script for using cluster machine
ST_num_frame	   = dic["num_frame"];         # number of frames in the first trajectory file
ST_num_atoms	   = dic["num_atoms"];         # number of atomes in the system
ST_num_files	   = dic["num_files"];         # number of files chosen
ST_num_ps	   = str(dic["num_ps"]);       # simulation time ps/frame
ST_trajectoryFile = [];                        # the list of trajectory files
ST_STATIC = ST_ANALYZER_HOME[0:len(ST_ANALYZER_HOME)-10];  # directory where ~/dJango_home/static

for i in range(len(dic["trajectory"])):
    ST_trajectoryFile.append(dic["trajectory"][i]);

# Identifying my parameters
myPara = get_myfunction(exe_file, dic);
fName = myPara[0];      # name of function except '.py'
pInfo = myPara[1];      # parameter information pInfo[0] = "number of parameters"
paras = myPara[2];      # actual parameters paras[0] contains 'the number of parameters'
para_pkey = myPara[3];  # primary key of parameter table contating this analzyer function. 
ST_rmodule   = "{0}{1}".format(exe_file[:len(exe_file)-3], ST_para_idx);  # running module name (e.g. box0)

#---------------------< assigned module specific parameters: fixed for every interface >---------------------------------
ST_num_paras = paras[0][0];		# pInfo[0] : number of parameters
ST_frmInt    = paras[1][ST_para_idx];	# pInfo[1] : Frame interval (list)
ST_outFile   = paras[2][ST_para_idx];	# pInfo[2] : output file name (list)

# Updating DB: running Job
# 0 - submit job
# 1 - Running job
# 2 - Error occurred
# 3 - Completed
stime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
conn = sqlite3.connect(ST_DB_FILE);
c    = conn.cursor();
for i in range(len(para_pkey)):
    query = """UPDATE gui_parameter SET status = "RUNNING" WHERE id = {}""".format(para_pkey[i]);
    c.execute(query);
    conn.commit();
conn.close();

# Update values into gui_outputs
conn = sqlite3.connect(ST_DB_FILE);
c    = conn.cursor();
# find gui_outputs related to current processing
jobName = "{}_{}".format(ST_rmodule, ST_outFile);
query = """SELECT id FROM gui_outputs WHERE job_id = {0} and name = "{1}" """.format(ST_job_pkey[0], jobName);
c.execute(query);
row = c.fetchone();
pk_output = row[0];     # primary key for gui_outputs
try:    
    query = """UPDATE gui_outputs SET status = "Running" WHERE id = {0}""".format(pk_output);
    c.execute(query);
    conn.commit();
    conn.close();
except:
    conn = sqlite3.connect(ST_DB_FILE);
    c    = conn.cursor();
    query = """UPDATE gui_outputs SET status = "Failed" WHERE id = {0}""".format(pk_output);
    c.execute(query);
    conn.commit();
    conn.close();


ST_out_dir = "{}{}".format(ST_OUTPUT_HOME, fName[1]);
# Create output directory
if not (os.path.isdir(ST_out_dir)):
    os.mkdir(ST_out_dir);


# -------- Writing input file for web-link
inFile = '{0}/input{1}.dat'.format(ST_out_dir, ST_para_idx);
fid_in = open(inFile, 'w');
strPara = "Name of Function: {}\n".format(exe_file);
strPara = strPara + "System information \n";
strPara = strPara + "\t- First trajectory contains {0} frames ({1} ps/frame)\n".format(ST_num_frame, ST_num_ps);
strPara = strPara + "\t- There are {0} trajectory file(s) and {1} atoms\n".format(ST_num_files, ST_num_atoms);
strPara = strPara + "Total {} parameters \n".format(int(paras[0][0])+3);
strPara = strPara + "\t- Base path: {}\n".format(ST_base_path);
strPara = strPara + "\t- Structure file: {}\n".format(ST_structure_file);
strPara = strPara + "\t- Trajectory files: \n"
tmp = "";
for trj in ST_trajectoryFile:
    tmp = tmp + "\t\t{}\n".format(trj);
    
strPara = strPara + tmp;

strPara = strPara + "\t- Job specific parameters: \n"
strPara = strPara + "\t\t{0}:{1}\n".format(pInfo[0], paras[0][0]);
tmp = "";
for i in range(len(pInfo)):
    if i > 0:
	tmp = tmp + "\t\t{0}:{1}\n".format(pInfo[i], paras[i][ST_para_idx]);
strPara = strPara + tmp;
strPara = strPara + "\nPBS: \n{}\n".format(ST_pbs);
fid_in.write(strPara);
fid_in.close();

"""
print '### Parameters ###'
print ST_para_idx 	 # parameter index for multiple jobs in a same form
print ST_para_pkeys   	 # gui_parameter primary keys obtained from job submitting
print ST_job_pkey        # gui_job primary keys obtained from job submitting
print ST_job_title       # gui_job title
print ST_prj_pkey        # gui_project primary key beloning to current job
print ST_prj_title       # gui_project title beloning to current job
print ST_SESSION_HOME    # default session home directory: ~/dJango_home/media/user_id
print ST_OUTPUT_HOME     # output directory given by user
print ST_PBS_HOME        # directory where PBS script is written (i.e. for using cluster)
print ST_ANALYZER_HOME   # directory where all analyzers are located (i.e. the location of this file)
print ST_MEDIA_HOME      # media directory ~/dJango_home/media
print ST_DB_FILE         # database file name and full location 
print ST_base_path       # base location of input files
print ST_path_output     # the location of output directory
print ST_path_python     # python path to run analyzers
print ST_structure_file  # full path of structure file (i.e. PDB, PSF, etc)
print ST_pbs             # PBS script for using cluster machine
print ST_num_frame	 # number of frames in the first trajectory file
print ST_num_atoms	 # number of atomes in the system
print ST_num_files	 # number of files chosen
print ST_num_ps	  	 # simulation time ps/frame
print ST_trajectoryFile  # the list of trajectory files
print ST_rmodule   	 # running module name (e.g. box0)
print ST_num_paras 	 # number of parameters
print ST_frmInt	  	 # Frame interval 
print ST_outFile   	 # output file name
"""
#---------------------< assigned global parameters >---------------------------------
# ST_para_idx 	  	# parameter index for multiple jobs in a same form
# ST_para_pkeys   	# gui_parameter primary keys obtained from job submitting
# ST_job_pkey        	# gui_job primary keys obtained from job submitting
# ST_job_title       	# gui_job title
# ST_prj_pkey        	# gui_project primary key beloning to current job
# ST_prj_title       	# gui_project title beloning to current job
# ST_SESSION_HOME    	# default session home directory: ~/dJango_home/media/user_id
# ST_OUTPUT_HOME     	# output directory given by user
# ST_PBS_HOME        	# directory where PBS script is written (i.e. for using cluster)
# ST_ANALYZER_HOME   	# directory where all analyzers are located (i.e. the location of this file)
# ST_STATIC		# static directory ~/dJango_home/static
# ST_MEDIA_HOME      	# media directory ~/dJango_home/media
# ST_DB_FILE         	# database file name and full location 
# ST_base_path       	# base location of input files
# ST_path_output     	# the location of output directory
# ST_path_python     	# python path to run analyzers
# ST_structure_file  	# full path of structure file (i.e. PDB, PSF, etc)
# ST_pbs             	# PBS script for using cluster machine
# ST_num_frame	  	# number of frames in the first trajectory file
# ST_num_atoms	  	# number of atomes in the system
# ST_num_files	  	# number of files chosen
# ST_num_ps	  	# simulation time ps/frame
# ST_trajectoryFile  	# the list of trajectory files
# ST_rmodule   	  	# running module name (e.g. box0)
# ST_num_paras 	  	# number of parameters
# ST_frmInt	  	# Frame interval 
# ST_outFile   	  	# output file name
# ST_out_dir		# output directory

###############################################################################################################################
######################################## PLEASE DO NOT MODIFY ABOVE THIST LINE!!!! ############################################
###############################################################################################################################

#---------------------< assigned module specific parameters: defined by users >---------------------------------
ST_num_atoms = paras[3][ST_para_idx];
ST_num_atoms = int(ST_num_atoms);
cntQry    = paras[4][ST_para_idx];
cntAxs    = paras[5][ST_para_idx];
selQry    = paras[6][ST_para_idx];

mk_single = int(paras[7][ST_para_idx]);

frmInt = stanalyzer.getSeqNumber(ST_frmInt);

#///////////////////////////////////////////////////////////////////////////
# Running actual job
#///////////////////////////////////////////////////////////////////////////
try:
    run = 1;
    if run:
	psf = '{0}{1}'.format(ST_base_path, ST_structure_file);
	timeStamp = [];         # time stamp for trajectory
	if mk_single == 1:
	    out_pdb = "{0}/{1}".format(ST_out_dir, ST_outFile);
	    fid_pdb = MDAnalysis.Writer(out_pdb, multiframe=True);
	    
	for idx in range(len(ST_trajectoryFile)):
	    
	    # turning on periodic boundary conditions
	    MDAnalysis.core.flags['use_periodic_selections'] = True
	    MDAnalysis.core.flags['use_KDTree_routines'] = False
	    
	    # reading trajectory
	    trj = '{0}{1}'.format(ST_base_path, ST_trajectoryFile[idx]);
	    u = Universe(psf, trj);
	    
	    if mk_single != 1:
		# define output name = dcd file name + user defined output file name
		out_pdb = "{0}/{1}_{2}.pdb".format(ST_out_dir, ST_trajectoryFile[idx], ST_outFile);
		fid_pdb = MDAnalysis.Writer(out_pdb, multiframe=True);
	    
	    # read based on frame
	    cnt_frm = 1;
	    
	    for ts in u.trajectory:
		#######################################################
		############## Write your code below ##################
		#######################################################
		if cnt_frm in frmInt:
		    #print "{} is in {}".format(cnt_frm, frmInt);
		    #======= Centeralization =========
		    if (cntQry != 'no') :
			stanalyzer.centerByRes3(ts, u, cntQry, 1, cntAxs);
		    #==================================
		    selAtoms = u.selectAtoms(selQry);
		    fid_pdb.write(selAtoms);
		
		cnt_frm = cnt_frm +1;
	    
	    if mk_single != 1:
		fid_pdb.close();
	
	if mk_single == 1:
	    fid_pdb.close();
	
	########################################
	###### Writing final output below ######
	########################################
	# pick one of output file names for DB usage
	out_file = out_pdb;
	
	########################################
	########### Drawing graphs #############
	########################################
	outImg  = '{0}{1}.png'.format(exe_file[:len(exe_file)-3], ST_para_idx);
	imgPath = "{0}/{1}".format(ST_out_dir, outImg);
	imgpdb = '{}/images/pdb.png'.format(ST_STATIC);
	subprocess.call(["cp", imgpdb, imgPath]);
	
	
	# gzip all reaults
	outZip = "{0}project_{1}_{2}{3}.tar.gz".format(ST_OUTPUT_HOME, ST_prj_pkey, fName[1], ST_para_idx);
	subprocess.call(["tar", "czf", outZip, ST_out_dir]);
	
	print "# gzip is okay"
	
	# Update values into gui_outputs
	conn = sqlite3.connect(ST_DB_FILE);
	c    = conn.cursor();
	query = """UPDATE gui_outputs SET status = "Complete", img="{0}", txt="{1}", gzip="{2}" WHERE id = {3}""".format(imgPath, out_file, outZip, pk_output);
	c.execute(query);
	conn.commit();
	conn.close();
	
    
    
    ######################################## PLEASE DO NOT MODIFY BELOW THIST LINE!!!! ############################################
    # update gui_parameter & gui_job table when job completed
    etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
    conn = sqlite3.connect(ST_DB_FILE);
    c    = conn.cursor();
    for i in range(len(para_pkey)):
	query = """UPDATE gui_parameter SET status = "COMPLETE" WHERE id = {0}""".format(para_pkey[i]);
	c.execute(query);
	conn.commit();
    
    # update gui_job if every status in gui_parameter are COMPLETE
    query = """SELECT DISTINCT(status) FROM gui_parameter WHERE job_id = {0}""".format(ST_job_pkey[0]);
    c.execute(query);
    ST = c.fetchall();
    
    if (len(ST) == 1) and (ST[0][0] == "COMPLETE"):
	etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
	query = """UPDATE gui_job SET status = "COMPLETE", etime = "{0}" WHERE id = {1}""".format(etime, ST_job_pkey[0]);
	c.execute(query);
	conn.commit();
    
	# making tar file
	outZip = "{0}project_{1}.tar.gz".format(ST_OUTPUT_HOME, ST_prj_pkey[0]);
	subprocess.call(["tar", "czf", outZip, ST_OUTPUT_HOME]);
    
    conn.close();

#///////////////////////////////////////////////////////////////////////////
# Finalizing  job
# -- Use following codes to make your own function
#///////////////////////////////////////////////////////////////////////////
except:
    # update gui_parameter & gui_job table when job failed
    print "ERROR OCCURRED!"
    etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
    conn = sqlite3.connect(ST_DB_FILE);
    c    = conn.cursor();
    for i in range(len(para_pkey)):
        query = """UPDATE gui_parameter SET status = "FAILED" WHERE id = {0}""".format(para_pkey[i]);
        c.execute(query);
        conn.commit();
    query = """UPDATE gui_job SET status = "INTERRUPTED" WHERE id = {0}""".format(ST_job_pkey[0]);
    c.execute(query);
    conn.commit();
    
    query = """UPDATE gui_outputs SET status = "Failed" WHERE id = {0}""".format(pk_output);
    c.execute(query);
    conn.commit();

    conn.close();

#///////////////////////////////////////////////////////
# print gui_job and gui_parameter table
#///////////////////////////////////////////////////////
stime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
conn = sqlite3.connect(ST_DB_FILE);
c    = conn.cursor();

#print "========= gui_job ==========="
query = "SELECT id, name, proj_id, anaz, status, output, stime, etime FROM gui_job";
c.execute(query);
job = c.fetchall();
query = "SELECT id, job_id, anaz, para, val, status FROM gui_parameter";
c.execute(query);
PR = c.fetchall();
conn.close();

