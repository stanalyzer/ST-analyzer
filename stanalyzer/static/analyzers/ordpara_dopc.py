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
num_atoms = paras[4][para_idx];			# pInfo[8] : Total number of atoms
num_atoms = int(num_atoms);

# defining the order of lipid tail
dnst_min = 2;
dnst_max = 16;
dnst_bin = 1;

# query selecting carbon tail
selQryC2 = "resname DPPE and name C2* and not name C21"; #SCD2
selQryC3 = "resname DXPE and name C3* and not name C31"; #SCD1
newQryC2 = "resname DPPE and (name C2* or name H*) and not name C21"; #SCD2
newQryC3 = "resname DXPE and (name C3* or name H*) and not name C31"; #SCD1

print "Okay I am in DPPE!!!!!";
print "AXIS = {}, {}".format(taxis, type(taxis));
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
	
	fid_out.write("# Range\tDensity\n");
	psf = '{0}{1}'.format(base_path, structure_file);
    
	cnt = 0;
	timeStamp = [];         # time stamp for trajectory
	
	# data based on trajectory
	DNST = [];
	STMP = [];
	for ibin in frange(dnst_min, dnst_max, dnst_bin):
	    DNST.append(0.0);
    
	# defining unit vector
	if taxis == 'X':
	    tx = 1.0;
	    ty = 0.0;
	    tz = 0.0;
	elif taxis == 'Y':
	    tx = 0.0;
	    ty = 1.0;
	    tz = 0.0;
	elif taxis == 'Z':
	    tx = 0.0;
	    ty = 0.0;
	    tz = 1.0;
	v2 = [tx, ty, tz];
	
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
		DegT = [];
		#tclock = cnt;
		cnt = cnt + 1;
		if (cnt % frmInt) == 0:
		    tmp_time = float(cnt) * float(num_ps) - float(num_ps);
		    STMP.append(tmp_time);
		    print "[{0}ps]selecting atoms...".format(tmp_time);
		    selAtoms = u.selectAtoms(selQryC3);
		    print "DONE!"
		    if len(selAtoms) > 1:
			# get coordinates
			CAs = selAtoms.names();
			ca0 = CAs[0:1];			# Ca index
			print "Ca0: {0}".format(ca0);
			
			c_ini = CAs[1];
			c_ini = int(c_ini[2:]);		# initial C index of tail
			
			c_end = CAs[len(selAtoms)-1];
			c_end = int(c_end[2:]);		# last C index of tail
    
			print "-- Order parameter: start:{0}, end:{1}".format(c_ini, c_end);
			newAtoms = u.selectAtoms(newQryC3);
			CRDs = newAtoms.coordinates();
    
			# get list of all atom names
			ANames = newAtoms.names();			
			for b_cnt in BIN:
			    tmp_Scd = 0.0;
			    num_Hs = 0.0;				# count number of Hyderogen at each carbon
			    c_cnt = int(b_cnt);
			    catom_idx = -1;
			    print "tmp_Scd={}".format(tmp_Scd);
			    
			    # defining current atoms 
			    Cxx = "{0}{1}".format(ca0[0], c_cnt);
			    HxR = "H{0}X".format(c_cnt);
			    HxS = "H{0}Y".format(c_cnt);
			    HxT = "H{0}Z".format(c_cnt);
			    print "c_cnt={0}, b_cnt={1}, Cxx={2}, HxR={3}, HxS={4}, HxT={5}".format(c_cnt, b_cnt, Cxx, HxR, HxS, HxT)
			    
			    # CA atom
			    if (Cxx in ANames):
				catom_idx = ANames.index(Cxx);
				CAx = CRDs[catom_idx][0];
				CAy = CRDs[catom_idx][1];
				CAz = CRDs[catom_idx][2];
				print "CAx={}, CAy={}, CAz={}".format(CAx, CAy, CAz);
				
			    # H-R atom
			    if (catom_idx >= 0) and (HxR in ANames):
				hratom_idx = ANames.index(HxR);
				HRx = CRDs[hratom_idx][0];
				HRy = CRDs[hratom_idx][1];
				HRz = CRDs[hratom_idx][2];
				
				x = HRx - CAx;
				y = HRy - CAy;
				z = HRz - CAz;
				v1 = [x, y, z];
				
				cosT = getCosT(v1, v2);
				tmp_Scd  += (3.0 * cosT * cosT - 1.0);
				num_Hs += 1.0;
			    
			    # H-S atom
			    if (catom_idx > 0) and (HxS in ANames):
				hsatom_idx = ANames.index(HxS);
				HSx = CRDs[hsatom_idx][0];
				HSy = CRDs[hsatom_idx][1];
				HSz = CRDs[hsatom_idx][2];
				
				x = HSx - CAx;
				y = HSy - CAy;
				z = HSz - CAz;
				
				v1 = [x, y, z];
    
				cosT = getCosT(v1, v2);
				tmp_Scd  += (3.0 * cosT * cosT - 1.0);
				num_Hs += 1.0;
    
			    # H-T atom
			    if  (catom_idx > 0) and (HxT in ANames):
				htatom_idx = ANames.index(HxT);
				HTx = CRDs[htatom_idx][0];
				HTy = CRDs[htatom_idx][1];
				HTz = CRDs[htatom_idx][2];
				
				x = HTx - CAx;
				y = HTy - CAy;
				z = HTz - CAz;
				v1 = [x, y, z];
    
				cosT = getCosT(v1, v2);
				tmp_Scd  += (3.0 * cosT * cosT - 1.0);
				num_Hs += 1.0;
			    
			    print "after tmp_Scd={}".format(tmp_Scd);
			    
			    # average Scd
			    total_Scd   = 0.5 * (tmp_Scd / num_Hs);
			    pos = bisect_left(BIN, b_cnt);
			    DNST[pos] += total_Scd;
			    print "+----------------------------------------------------------+";
    
    
    
	# Write down results
	finalDNST = [];
	for i in DNST:
	    tmp = abs(i/len(STMP));	# average through trajectory
	    finalDNST.append(tmp);
	
	# Writing final output
	for i in range(len(finalDNST)):
	    outStr = "{0}\t{1}\n".format(BIN[i], finalDNST[i]);
	    fid_out.write(outStr);
	fid_out.close()
    
	# -------- Drawing graphs
	# Writing Gnuplot script
	outScr = '{0}/gplot{1}.p'.format(out_dir, para_idx);
	outImg  = '{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
	imgPath = "{0}/{1}".format(out_dir, outImg);
	fid_out = open(outScr, 'w');
	gScript = "set terminal png\n";
	gScript = gScript + "set xlabel 'Carbon Index'\n";
	gScript = gScript + "set ylabel 'S_CD'\n";
	gScript = gScript + "set output '{0}'\n".format(imgPath);
	gScript = gScript + """plot "{0}/{1}" using 1:2 title "DOPC" with lines lw 3\n""".format(out_dir, outFile);
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

