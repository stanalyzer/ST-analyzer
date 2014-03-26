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
ST_frmInt    = int(ST_frmInt);
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
print '### Parameters in RDF1D.py ###'
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
cntQry	= paras[4][ST_para_idx];		# bin size
cntAxs	= paras[5][ST_para_idx];		# bin size
binsize = paras[6][ST_para_idx];		# bin size
binsize = int(binsize);
selQry1 = paras[7][ST_para_idx];		# selection query 1
selQry2 = paras[8][ST_para_idx];		# selection query 2

sysX    = paras[9][ST_para_idx];		# system size X
sysX 	= round(float(sysX));

sysY    = paras[10][ST_para_idx];		# system size Y
sysY 	= round(float(sysY));

sysZ    = paras[11][ST_para_idx];		# system size Z
sysZ 	= round(float(sysZ));

augX = sysX + 10.0;
augY = sysY + 10.0;
augZ = sysZ + 10.0;

# defining top and bottom leaflet 
if cntAxs == 'x':
    top_cut = sysX * 0.165;
    btm_cut = -1.0 * top_cut;
    top_selQry1 = "{} and (prop x > {})".format(selQry1, top_cut);
    top_selQry2 = "{} and (prop x > {})".format(selQry2, top_cut);
    btm_selQry1 = "{} and (prop x < {})".format(selQry1, btm_cut);
    btm_selQry2 = "{} and (prop x < {})".format(selQry2, btm_cut);
elif cntAxs == 'y':
    top_cut = sysY * 0.165;
    btm_cut = -1.0 * top_cut;
    top_selQry1 = "{} and (prop y > {})".format(selQry1, top_cut);
    top_selQry2 = "{} and (prop y > {})".format(selQry2, top_cut);
    btm_selQry1 = "{} and (prop y < {})".format(selQry1, btm_cut);
    btm_selQry2 = "{} and (prop y < {})".format(selQry2, btm_cut);
else:
    top_cut = sysZ * 0.165;
    btm_cut = -1.0 * top_cut;
    top_selQry1 = "{} and (prop z > {})".format(selQry1, top_cut);
    top_selQry2 = "{} and (prop z > {})".format(selQry2, top_cut);
    btm_selQry1 = "{} and (prop z < {})".format(selQry1, btm_cut);
    btm_selQry2 = "{} and (prop z < {})".format(selQry2, btm_cut);

#----- creating BIN for distance density
min_bin = 1.0;
max_bin = math.sqrt(augX * augX + augY * augY + augY * augY);

BIN = [];
for i in stanalyzer.frange(min_bin, max_bin, binsize):
    BIN.append(i);

#----- calculating denominator for density distribution
DEN = []; topDNST = []; btmDNST = []; 
pre_out = 0.0;
for i in range(len(BIN)):
    if i == 0:
	sp_out = 4.0/3.0 * math.pi * (BIN[i]**3);
	sp_in  = 0.0;
	vol = sp_out - sp_in;
	pre_out = sp_out;
    else:
	sp_out = 4.0/3.0 * math.pi * (BIN[i]**3);
	sp_in = pre_out;
	vol = sp_out - sp_in;
	pre_out = sp_out;
    DEN.append(vol);
    topDNST.append(0.0);
    btmDNST.append(0.0);

print "Okay running...."

#///////////////////////////////////////////////////////////////////////////
# Running actual job
#///////////////////////////////////////////////////////////////////////////
try:
    run = 1;
    if run:
	psf = '{0}{1}'.format(ST_base_path, ST_structure_file);

	for idx in range(len(ST_trajectoryFile)):
	    # turning on periodic boundary conditions
	    MDAnalysis.core.flags['use_periodic_selections'] = True
	    MDAnalysis.core.flags['use_KDTree_routines'] = False
	    
	    # reading trajectory
	    trj = '{0}{1}'.format(ST_base_path, ST_trajectoryFile[idx]);
	    u = Universe(psf, trj);

	    # calculating dimensions
	    sysD = u.dimensions;
	    csize_x = sysD[0];	# current system size of X axis
	    csize_y = sysD[1];
	    csize_z = sysD[2];

	    allAtoms = u.selectAtoms('all');
	    
	    cnt = 0;
	    num_frm = 0.0;
	    
	    for ts in u.trajectory:
		# flag for top and bottom leaflet
		top_flag = 0; btm_flag = 0;
		
		# initializing array
		tmp_topDNST = [0.0] * len(BIN);
		tmp_btmDNST = [0.0] * len(BIN);

		cnt = cnt + 1;
		
		if (cnt % ST_frmInt) == 0:
		    num_frm = num_frm + 1.0;
		    
		    #######################################################
		    ############## Write your code below ##################
		    #######################################################
		    
		    #======= Centeralization =========
		    if (cntQry != 'no') :
			stanalyzer.centerByRes2(ts, u, cntQry, 1, cntAxs); # 1st residue is always chosen for centering membrane
		    #==================================
    
		    #############################################
		    #### Calculating top and bottom layer #######
		    #############################################
		    # get atom coordinates of each selection query
			    
		    ####################
		    #----> top leaflet
		    ####################
		    #---> Source Atoms
		    top_sAtoms = u.selectAtoms(top_selQry1);
		    top_sRes = top_sAtoms.resnames();
		    top_sIDs = top_sAtoms.resids();
		    top_sQry = '';
		    for i in range(len(top_sRes)):
			if i == 0:
			    top_sQry = '((resname {} and resid {})'.format(top_sRes[i], top_sIDs[i]);
			elif (i > 0) and  (i < len(top_sRes)-1):
			    top_sQry = '{} or (resname {} and resid {})'.format(top_sQry, top_sRes[i], top_sIDs[i]);
			else:
			    top_sQry = '{} or (resname {} and resid {})) and ({})'.format(top_sQry, top_sRes[i], top_sIDs[i], selQry1);
		    #print "##### top_sAtoms: {}".format(top_sQry);
		    top_sAtoms = u.selectAtoms(top_sQry);
		    if len(top_sAtoms.resnames()) > 0:
			top_sCRDs = top_sAtoms.coordinates();
			top_flag = top_flag + 1;
		    
		    #---> Target Atoms
		    top_tAtoms = u.selectAtoms(top_selQry2);
		    top_tRes = top_tAtoms.resnames();
		    top_tIDs = top_tAtoms.resids();
		    top_tQry = '';
		    for i in range(len(top_tRes)):
			if i == 0:
			    top_tQry = '((resname {} and resid {})'.format(top_tRes[i], top_tIDs[i]);
			elif (i > 0) and  (i < len(top_tRes)-1):
			    top_tQry = '{} or (resname {} and resid {})'.format(top_tQry, top_tRes[i], top_tIDs[i]);
			else:
			    top_tQry = '{} or (resname {} and resid {})) and ({})'.format(top_tQry, top_tRes[i], top_tIDs[i], selQry2);
		    #print "##### top_tAtoms: {}".format(top_tQry);
		    top_tAtoms = u.selectAtoms(top_tQry);
		    if len(top_tAtoms.resnames()) > 0:
			top_tCRDs = top_tAtoms.coordinates();
			top_flag = top_flag + 1;
		    
		    ######################
		    #----> bottom leaflet
		    ######################
		    #---> Source Atoms
		    btm_sAtoms = u.selectAtoms(btm_selQry1);
		    btm_sRes = btm_sAtoms.resnames();
		    btm_sIDs = btm_sAtoms.resids();
		    btm_sQry = '';
		    for i in range(len(btm_sRes)):
			if i == 0:
			    btm_sQry = '((resname {} and resid {})'.format(btm_sRes[i], btm_sIDs[i]);
			elif (i > 0) and  (i < len(btm_sRes)-1):
			    btm_sQry = '{} or (resname {} and resid {})'.format(btm_sQry, btm_sRes[i], btm_sIDs[i]);
			else:
			    btm_sQry = '{} or (resname {} and resid {})) and ({})'.format(btm_sQry, btm_sRes[i], btm_sIDs[i], selQry1);
		    #print "##### btm_sAtoms: {}".format(btm_sQry);
		    btm_sAtoms = u.selectAtoms(btm_sQry);
		    if len(btm_sAtoms.resnames()) > 0:
			btm_sCRDs = btm_sAtoms.coordinates();
			btm_flag = btm_flag + 1;
		    
		    #---> Target Atoms
		    btm_tAtoms = u.selectAtoms(btm_selQry2);
		    btm_tRes = btm_tAtoms.resnames();
		    btm_tIDs = btm_tAtoms.resids();
		    btm_tQry = '';
		    for i in range(len(btm_tRes)):
			if i == 0:
			    btm_tQry = '((resname {} and resid {})'.format(btm_tRes[i], btm_tIDs[i]);
			elif (i > 0) and  (i < len(btm_tRes)-1):
			    btm_tQry = '{} or (resname {} and resid {})'.format(btm_tQry, btm_tRes[i], btm_tIDs[i]);
			else:
			    btm_tQry = '{} or (resname {} and resid {})) and ({})'.format(btm_tQry, btm_tRes[i], btm_tIDs[i], selQry2);
		    #print "##### btm_tAtoms: {}".format(btm_tQry);
		    btm_tAtoms = u.selectAtoms(btm_tQry);
		    if len(btm_tAtoms.resnames()) > 0:
			btm_tCRDs  = btm_tAtoms.coordinates();
			btm_flag = btm_flag + 1;
		    
		    # for all atoms to get the coordinates
		    allAtoms   = u.selectAtoms('all');
		    allCRDs	   = allAtoms.coordinates();
		    
		    # considering periodic boundary for distance calculation
		    a_mid_x = 0.5 * (max(allCRDs[:,0]) + min(allCRDs[:,0]));
		    a_mid_y = 0.5 * (max(allCRDs[:,1]) + min(allCRDs[:,1]));
		    a_mid_z = 0.5 * (max(allCRDs[:,2]) + min(allCRDs[:,2]));
		    
    
		    # calculating distance & density calculation for top layer
		    # we calculated each of them due to uneven bilayer
		    #print "##### top_flag = {}#####".format(top_flag);
		    #print "##### btm_flag = {}#####".format(btm_flag);
		    
		    if top_flag == 2:
			#print "---> top leaflet calculation"
			for i in range(len(top_sCRDs)):
			    for j in range(len(top_tCRDs)):
				# adjusting X axis
				if (top_tCRDs[j,0] >= a_mid_x):
				    adjX = top_tCRDs[j,0] - csize_x;
				else:
				    adjX = top_tCRDs[j,0] + csize_x;
				    
				# adjusting Y axis
				if (top_tCRDs[j,1] >= a_mid_y):
				    adjY = top_tCRDs[j,1] - csize_y;
				else:
				    adjY = top_tCRDs[j,1] + csize_y;
				    
				# adjusting Z axis
				if (top_tCRDs[j,2] >= a_mid_z):
				    adjZ = top_tCRDs[j,2] - csize_z;
				else:
				    adjZ = top_tCRDs[j,2] + csize_z;
				
				tmpCRDs = np.array([adjX, adjY, adjZ]);
				
				dist1 = np.sqrt(np.sum((top_sCRDs[i] - top_tCRDs[j]) ** 2));
				dist2 = np.sqrt(np.sum((top_sCRDs[i] - tmpCRDs) ** 2));
				dist  = min([dist1, dist2]);
				
				tmp = np.array(stanalyzer.count_intervals2(dist, BIN));
				tmp_topDNST = tmp_topDNST + tmp;
    
		    if btm_flag == 2:
			#print "---> bottom leaflet calculation"
			for i in range(len(btm_sCRDs)):
			    for j in range(len(btm_tCRDs)):
				# adjusting X axis
				if (btm_tCRDs[j,0] >= a_mid_x):
				    adjX = btm_tCRDs[j,0] - csize_x;
				else:
				    adjX = btm_tCRDs[j,0] + csize_x;
				    
				# adjusting Y axis
				if (btm_tCRDs[j,1] >= a_mid_y):
				    adjY = btm_tCRDs[j,1] - csize_y;
				else:
				    adjY = btm_tCRDs[j,1] + csize_y;
				    
				# adjusting Z axis
				if (btm_tCRDs[j,2] >= a_mid_z):
				    adjZ = btm_tCRDs[j,2] - csize_z;
				else:
				    adjZ = btm_tCRDs[j,2] + csize_z;
				
				tmpCRDs = np.array([adjX, adjY, adjZ]);
				
				dist1 = np.sqrt(np.sum((btm_sCRDs[i] - btm_tCRDs[j]) ** 2));
				dist2 = np.sqrt(np.sum((btm_sCRDs[i] - tmpCRDs) ** 2));
				dist  = min([dist1, dist2]);
				
				tmp = np.array(stanalyzer.count_intervals2(dist, BIN));
				tmp_btmDNST = tmp_btmDNST + tmp;
	    
	    # to prevent overflow tmp_btmDNST is used
	    if top_flag == 2:
		topDNST = topDNST + tmp_topDNST / num_frm;
	    
	    if btm_flag == 2:
		btmDNST = btmDNST + tmp_btmDNST / num_frm;
	    
	if top_flag == 2:
	    topDNST = topDNST / len(ST_trajectoryFile);
	    topDNST = topDNST / DEN;
	
	if btm_flag == 2:
	    btmDNST = btmDNST / len(ST_trajectoryFile);
	    btmDNST = btmDNST / DEN;
			
	########################################
	###### Writing final output below ######
	########################################
	ftop = '{0}/top_{1}'.format(ST_out_dir, ST_outFile);
	fid_ftop = open(ftop, 'w');
	fbtm = '{0}/btm_{1}'.format(ST_out_dir, ST_outFile);
	fid_fbtm = open(fbtm, 'w');
    
	# pick one of output file names for DB usage
	out_file = ftop;
	
	# writing comments '#' is considered as a comment in GNUPlot	
	cmt = "#BIN\tDistribution\n";
	fid_ftop.write(cmt);
	fid_fbtm.write(cmt);
	for i in range(len(topDNST)):
	    out_data_top = "{}\t{}\n".format(BIN[i], topDNST[i]);
	    out_data_btm = "{}\t{}\n".format(BIN[i], btmDNST[i]);
	    fid_ftop.write(out_data_top);
	    fid_fbtm.write(out_data_btm);
	fid_ftop.close();
	fid_fbtm.close();
	
	########################################
	########### Drawing graphs #############
	########################################
	# Writing Gnuplot script

	outScr = '{0}/gplot{1}.gpl'.format(ST_out_dir, ST_para_idx);
	outImg  = '{0}{1}.png'.format(exe_file[:len(exe_file)-3], ST_para_idx);
	imgPath = "{0}/{1}".format(ST_out_dir, outImg);
	fid_out = open(outScr, 'w');
	gScript = """set terminal png enhanced \n""";
	gScript = gScript + "set output '{0}'\n".format(imgPath);
	gScript = gScript + "set multiplot layout 2, 1 title 'RDF 1D'\n";
	
	gScript = gScript + """set xtics offset 0,0.5; unset xlabel\n""";
	gScript = gScript + """set label 1 'Top' at graph 0.01, 0.95 font ',8'\n""";
	gScript = gScript + """set ylabel 'Density' offset 0,-8\n""";
	gScript = gScript + """plot "{0}" using 1:2 title "RDF 1D" with lines lw 3\n""".format(ftop);
	
	gScript = gScript + """set xtics offset 0,0.5; unset xlabel; unset ylabel\n""";
	gScript = gScript + """set label 1 'Bottom' at graph 0.01, 0.95 font ',8'\n""";
	gScript = gScript + """set xlabel 'Distance' offset 0,1\n""";
	gScript = gScript + """plot "{0}" using 1:2 title "RDF 1D" with lines lw 3\n""".format(fbtm);
	
	gScript = gScript + "unset multiplot\n";
	fid_out.write(gScript);
	fid_out.close();
	
	# Drawing graph with gnuplot
	subprocess.call(["gnuplot", outScr]);
	
	# gzip all reaults
	outZip = "{0}project_{1}_{2}{3}.tar.gz".format(ST_OUTPUT_HOME, ST_prj_pkey, fName[1], ST_para_idx);
	subprocess.call(["tar", "czf", outZip, ST_out_dir]);
	
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

