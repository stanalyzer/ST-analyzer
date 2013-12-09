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
para_idx = int(sys.argv[2]);         # parameter index for multiple jobs in a same form
#out_file = sys.argv[2];

print '### execution file: {} ###'.format(exe_file);
print '\t- input file: {}'.format(in_file);
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
pdb_file  	= dic["pdb_file"];    	    # full path of structure file (i.e. PDB)
pbs             = dic["pbs"];               # PBS script for using cluster machine
num_frame	= dic["num_frame"];         # number of frames in the first trajectory file
num_atoms	= dic["num_atoms"];         # number of atomes in the system
num_files	= dic["num_files"];         # number of files chosen
num_ps	        = str(dic["num_ps"]);       # simulation time ps/frame
trajectoryFile = [];                        # the list of trajectory files
for i in range(len(dic["trajectory"])):
    trajectoryFile.append(dic["trajectory"][i]);

# Identifying my parameters
myPara = get_myfunction(exe_file, dic);
fName = myPara[0];      # name of function except '.py'
pInfo = myPara[1];      # parameter information pInfo[0] = "number of parameters"
paras = myPara[2];      # actual parameters paras[0] contains 'the number of parameters'
para_pkey = myPara[3];  # primary key of parameter table contating this analzyer function. 
rmodule   = "{0}{1}".format(exe_file[:len(exe_file)-3], para_idx);  # running module name (e.g. box0)

#---------------------< assigned module specific parameters: fixed for every interface >---------------------------------
num_paras = paras[0][0];			# pInfo[0] : number of parameters
frmInt	  = paras[1][para_idx];			# pInfo[1] : Frame interval (list)
frmInt	  = int(frmInt);
outFile   = paras[2][para_idx];			# pInfo[2] : output file name (list)

#print "NAME OF FUNCTION: {}".format(fName);
#print "PARAMETER INFO: {}".format(pInfo);
#print "PARMETERS: {}".format(paras);
#print "PARAMETER PRIMARY KEYS: {}".format(para_pkey);

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

# Update values into gui_outputs
conn = sqlite3.connect(DB_FILE);
c    = conn.cursor();
# find gui_outputs related to current processing
jobName = "{}_{}".format(rmodule, outFile);
query = """SELECT id FROM gui_outputs WHERE job_id = {0} and name = "{1}" """.format(job_pkey[0], jobName);
c.execute(query);
row = c.fetchone();
print query
pk_output = row[0];     # primary key for gui_outputs
try:    
    query = """UPDATE gui_outputs SET status = "Running" WHERE id = {0}""".format(pk_output);
    c.execute(query);
    conn.commit();
    conn.close();
except:
    conn = sqlite3.connect(DB_FILE);
    c    = conn.cursor();
    query = """UPDATE gui_outputs SET status = "Failed" WHERE id = {0}""".format(pk_output);
    c.execute(query);
    conn.commit();
    conn.close();


out_dir = "{}{}".format(OUTPUT_HOME, fName[1]);
# Create output directory
if not (os.path.isdir(out_dir)):
    #print "Creating directory into {}".format(out_dir)
    os.mkdir(out_dir);


# -------- Writing input file for web-link
#print "list of PARAMETERS: "
inFile = '{0}/input{1}.dat'.format(out_dir, para_idx);
fid_in = open(inFile, 'w');
strPara = "Name of Function: {}\n".format(exe_file);
strPara = strPara + "System information \n";
strPara = strPara + "\t- First trajectory contains {0} frames ({1} ps/frame)\n".format(num_frame, num_ps);
strPara = strPara + "\t- There are {0} trajectory file(s) and {1} atoms\n".format(num_files, num_atoms);
strPara = strPara + "Total {} parameters \n".format(int(paras[0][0])+3);
strPara = strPara + "\t- Base path: {}\n".format(base_path);
strPara = strPara + "\t- Structure file: {}\n".format(structure_file);
strPara = strPara + "\t- Trajectory files: \n"
tmp = "";
for trj in trajectoryFile:
    tmp = tmp + "\t\t{}\n".format(trj);
    
strPara = strPara + tmp;

strPara = strPara + "\t- Job specific parameters: \n"
strPara = strPara + "\t\t{0}:{1}\n".format(pInfo[0], paras[0][0]);
tmp = "";
for i in range(len(pInfo)):
    if i > 0:
	tmp = tmp + "\t\t{0}:{1}\n".format(pInfo[i], paras[i][para_idx]);
strPara = strPara + tmp;
strPara = strPara + "\nPBS: \n{}\n".format(pbs);
fid_in.write(strPara);
fid_in.close();

#print strPara

#---------------------< assigned global parameters >---------------------------------
# para_idx 	  = int(sys.argv[2]);         # parameter index for multiple jobs in a same form
# para_pkeys      = dic["para_pkeys"];        # gui_parameter primary keys obtained from job submitting
# job_pkey        = dic["job_pkey"];          # gui_job primary keys obtained from job submitting
# job_title       = dic["job_title"];         # gui_job title
# prj_pkey        = dic["pkey"];              # gui_project primary key beloning to current job
# prj_title       = dic["ptitle"];            # gui_project title beloning to current job
# SESSION_HOME    = dic["session_home"];      # default session home directory: ~/dJango_home/media/user_id
# OUTPUT_HOME     = dic["output_home"];       # output directory given by user
# PBS_HOME        = dic["pbs_home"];          # directory where PBS script is written (i.e. for using cluster)
# ANALYZER_HOME   = dic["analyzer_home"];     # directory where all analyzers are located (i.e. the location of this file)
# MEDIA_HOME      = dic["media_home"];        # media directory ~/dJango_home/media
# DB_FILE         = dic["dbName"];            # database file name and full location 
# base_path       = dic["base_path"];         # base location of input files
# path_output     = dic["path_output"];       # the location of output directory
# path_python     = dic["path_python"];       # python path to run analyzers
# structure_file  = dic["structure_file"];    # full path of structure file (i.e. PDB, PSF, etc)
# pbs             = dic["pbs"];               # PBS script for using cluster machine
# num_frame	  = dic["num_frame"];         # number of frames in the first trajectory file
# num_atoms	  = dic["num_atoms"];         # number of atomes in the system
# num_files	  = dic["num_files"];         # number of files chosen
# num_ps	  = str(dic["num_ps"]);       # simulation time ps/frame
# trajectoryFile  = [];                        # the list of trajectory files
# rmodule   	  = "{0}{1}".format(exe_file[:len(exe_file)-3], para_idx);  # running module name (e.g. box0)


###############################################################################################################################
######################################## PLEASE DO NOT MODIFY ABOVE THIST LINE!!!! ############################################
###############################################################################################################################

#---------------------< assigned module specific parameters: defined by users >---------------------------------
taxis	  = paras[3][para_idx];			# pInfo[3] : target axis 'ALL' 'X' 'Y' 'Z'
dnst_min  = paras[4][para_idx];			# pInfo[4] : Density min
dnst_min = float(dnst_min);

dnst_max  = paras[5][para_idx];			# pInfo[5] : Density max
dnst_max = float(dnst_max);

dnst_bin  = paras[6][para_idx];			# pInfo[6] : Bin size
dnst_bin = float(dnst_bin);

selQry = paras[7][para_idx];			# pInfo[7] : Query

num_atoms = paras[8][para_idx];			# pInfo[8] : Total number of atoms
num_atoms = int(num_atoms);

cntQry = paras[9][para_idx];			# pInfo[9] : Centering Query
cntAxs = paras[10][para_idx];			# pInfo[10]: Centering axis

print "AXIS = {}, {}".format(taxis, type(taxis));
print "MIN = {}, {}".format(dnst_min, type(dnst_min));
print "MAX = {}, {}".format(dnst_max, type(dnst_max));
print "BIN = {}, {}".format(dnst_bin, type(dnst_bin));
print "QUERY = {}, {}".format(selQry, type(selQry));
print "Centering QUERY = {}, {}".format(cntQry, type(cntQry));
print "Centering Axis = {}, {}".format(cntAxs, type(cntAxs));
print "Total # atoms = {}, {}".format(num_atoms, type(num_atoms));

#dummy = raw_input("Pause: ");

#///////////////////////////////////////////////////////////////////////////
# Running actual job
#///////////////////////////////////////////////////////////////////////////

try:
    run = 1;
    if run:
	out_file = '{0}/{1}'.format(out_dir, outFile);
	fid_out = open(out_file, 'w')
	cmt = "# Min={0}, Max={1}, Bin size = {2}\n".format(dnst_min, dnst_max, dnst_bin);
	fid_out.write(cmt)
	
	# describing BIN
	cmt = "# Bin ranges\n#"
	fid_out.write(cmt);
	
	# calculating BIN range
	BIN = [];
	for ibin in frange(dnst_min, dnst_max, dnst_bin):
	    BIN.append(ibin);
	    cmt = "\t{}".format(ibin);
	    fid_out.write(cmt);
	fid_out.write("\n");
	
	fid_out.write("# Range\tnumber density\tmass density\n");
	psf = '{0}{1}'.format(base_path, structure_file);
    
	cnt = 0;
	timeStamp = [];         # time stamp for trajectory
	
	# data based on trajectory
	DNST = [];
	mDNST = [];  # mass DNST
	STMP = [];
	
	for ibin in frange(dnst_min, dnst_max, dnst_bin):
	    DNST.append(0.0);
	    mDNST.append(0.0);
	    
	for idx in range(len(trajectoryFile)):
	    
	    # turning on periodic boundary conditions
	    MDAnalysis.core.flags['use_periodic_selections'] = True
	    MDAnalysis.core.flags['use_KDTree_routines'] = False
	    
	    # reading trajectory
	    trj = '{0}{1}'.format(base_path, trajectoryFile[idx]);
	    print 'Reading PSF: ' + psf
	    print 'Reading DCD: ' + trj
	    u = Universe(psf, trj);
	    print '{0} is done!'.format(idx);
	    
	    # read based on frame
	    for ts in u.trajectory:
		#tclock = cnt;
		cnt = cnt + 1;
		if (cnt % frmInt) == 0:
		    #======= Centeralization =========
		    if (cntQry != 'no') :
			#stanalyzer.centerByCOM(ts, u, cntQry);
			stanalyzer.centerByRes(ts, u, cntQry, 1, cntAxs); # 1st residue is always chosen for centering membrane
		    #==================================
		    
		    tmp_time = float(cnt) * float(num_ps) - float(num_ps);
		    STMP.append(tmp_time);
		    print "[{0}ps]selecting atoms...".format(tmp_time);
		    selAtoms = u.selectAtoms(selQry);
		    print "DONE!"
		    if len(selAtoms) > 1:
			# get coordinates
			CRDs = selAtoms.coordinates();
			MASS = selAtoms.masses();
			
			# get coordinate X, Y, Z
			crdX = []; crdY = []; crdZ =[];
			for iatom in range(len(CRDs)):
			    crdX.append(CRDs[iatom][0]);
			    crdY.append(CRDs[iatom][1]);
			    crdZ.append(CRDs[iatom][2]);
			
			sizeX = max(crdX) - min(crdX);
			sizeY = max(crdY) - min(crdY);
			sizeZ = max(crdZ) - min(crdZ);
			
			if taxis == 'X':
			    tcrd = crdX;
			    bin_vol = sizeY * sizeZ * dnst_bin;
			    print 'crdX is selected!';
			elif taxis == 'Y':
			    tcrd = crdY;
			    bin_vol = sizeX * sizeZ * dnst_bin;
			    print 'crdY is selected!';
			else:
			    tcrd = crdZ;
			    bin_vol = sizeX * sizeY * dnst_bin;
			    print 'crdZ is selected!';
			    
			# counting items based on BIN interval
			tmpDNST = count_intervals(tcrd, BIN);
			tmpDNST2 = count_intervals_mass(tcrd, MASS, BIN);
			
			# sort key based on numbering*
			intKEY = [];
			for key in sorted(tmpDNST.iterkeys()):
			    intKEY.append(float(key));
			intKEY.sort();

			intKEY2 = [];
			for key in sorted(tmpDNST2.iterkeys()):
			    intKEY2.append(float(key));
			intKEY2.sort();
			
			# reordering dictionary based on sorted key
			for key in intKEY:
			    idx_bin = BIN.index(key);
			    DNST[idx_bin] += float(tmpDNST[key]) / float(bin_vol); # normalizing based on volume
			    
			for key in intKEY2:
			    idx_bin = BIN.index(key);
			    mDNST[idx_bin] += float(tmpDNST2[key]) / float(bin_vol); # normalizing based on volume

	# Write down results
	finalDNST = [];
	for i in DNST:
	    tmp = i/len(STMP);
	    finalDNST.append(tmp);

	finalDNST2 = [];
	for i in mDNST:
	    tmp = i/len(STMP);
	    finalDNST2.append(tmp);
	
	# Writing final output
	for i in range(len(finalDNST)):
	    outStr = "{0}\t{1}\t{2}\n".format(BIN[i], finalDNST[i], finalDNST2[i]);
	    fid_out.write(outStr);
	fid_out.close()
    
	# -------- Drawing graphs
	# Writing Gnuplot script
	outScr = '{0}/gplot{1}.p'.format(out_dir, para_idx);
	outImg  = '{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
	imgPath = "{0}/{1}".format(out_dir, outImg);
	fid_out = open(outScr, 'w');
	gScript = "set terminal png\n";
	gScript = gScript + "set xlabel 'range'\n";
	gScript = gScript + "set ylabel 'density'\n";
	gScript = gScript + "set output '{0}'\n".format(imgPath);
	gScript = gScript + """plot "{0}/{1}" using 1:3 title "mass density" with lines lw 3\n""".format(out_dir, outFile);
	fid_out.write(gScript);
	fid_out.close()
	
	# Drawing graph with gnuplot
	subprocess.call(["gnuplot", outScr]);
	
	# gzip all reaults
	outZip = "{0}project_{1}_{2}{3}.tar.gz".format(OUTPUT_HOME, prj_pkey, fName[1], para_idx);
	subprocess.call(["tar", "czf", outZip, out_dir]);
    
	# Update values into gui_outputs
	conn = sqlite3.connect(DB_FILE);
	c    = conn.cursor();
	query = """UPDATE gui_outputs SET status = "Complete", img="{0}", txt="{1}", gzip="{2}" WHERE id = {3}""".format(imgPath, out_file, outZip, pk_output);
	c.execute(query);
	conn.commit();
	conn.close();
	#print query
    
    
    
    
    ######################################## PLEASE DO NOT MODIFY BELOW THIST LINE!!!! ############################################
    # update gui_parameter & gui_job table when job completed
    etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
    conn = sqlite3.connect(DB_FILE);
    c    = conn.cursor();
    for i in range(len(para_pkey)):
	query = """UPDATE gui_parameter SET status = "COMPLETE" WHERE id = {0}""".format(para_pkey[i]);
	#print query
	c.execute(query);
	conn.commit();
    
    # update gui_job if every status in gui_parameter are COMPLETE
    query = """SELECT DISTINCT(status) FROM gui_parameter WHERE job_id = {0}""".format(job_pkey[0]);
    c.execute(query);
    ST = c.fetchall();
    
    #print query;
    #print "number status = {}".format(len(ST));
    #for item in ST:
    #    print "{0}".format(item[0]);
    
    
    if (len(ST) == 1) and (ST[0][0] == "COMPLETE"):
	etime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
	query = """UPDATE gui_job SET status = "COMPLETE", etime = "{0}" WHERE id = {1}""".format(etime, job_pkey[0]);
	c.execute(query);
	conn.commit();
    
	# making tar file
	outZip = "{0}project_{1}.tar.gz".format(OUTPUT_HOME, prj_pkey[0]);
	subprocess.call(["tar", "czf", outZip, OUTPUT_HOME]);
    
	# Inserting compressed tar file for all submitted jobs
	#final_title = "[** All JOBs **] {0}".format(job_title);
	#query = """INSERT INTO gui_outputs (job_id, name, img, txt, gzip) VALUES ({0}, "{1}", "{2}", "{3}", "{4}")""".format(job_pkey[0], final_title, '', '', outZip);
	#c.execute(query);
	#conn.commit();
	
    conn.close();


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
    
    query = """UPDATE gui_outputs SET status = "Failed" WHERE id = {0}""".format(pk_output);
    c.execute(query);
    conn.commit();

    conn.close();

#///////////////////////////////////////////////////////
# print gui_job and gui_parameter table
#///////////////////////////////////////////////////////
stime = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
conn = sqlite3.connect(DB_FILE);
c    = conn.cursor();

#print "========= gui_job ==========="
query = "SELECT id, name, proj_id, anaz, status, output, stime, etime FROM gui_job";
#print "ID\tTITLE\tPROJ_ID\tANALYZER\tSTATUS\tOUTPUT\tSTART\tEND";
c.execute(query);
job = c.fetchall();
#print "Final idx= {}".format(job[len(job)-1][0]);
#for item in job:
#    print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(item[0], item[1], item[2], item[3], item[4], item[5], item[6], item[7]);

#print "========= gui_parameter ==========="
query = "SELECT id, job_id, anaz, para, val, status FROM gui_parameter";
#print "ID\tJOB_ID\tANALYZER\tPARAMETER\tVALUE\tSTATUS";
c.execute(query);
PR = c.fetchall();
#print "Final idx= {}".format(PR[len(PR)-1][0]);
#for item in PR:
#    print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(item[0], item[1], item[2], item[3], item[4], item[5]);

conn.close();

