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
#import lipidArea

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
print paras
segID 	= paras[3][para_idx];			# pInfo[3] : Total number of atoms
resName	= paras[4][para_idx];
selQry  = "segid " + segID + " and resname " + resName + " and (" + paras[5][para_idx] + ")";

print "Okay I am in sterol_tilt_ring.py!!!!!";
print "selQry = {}".format(selQry);

#///////////////////////////////////////////////////////////////////////////
# Running actual job
#///////////////////////////////////////////////////////////////////////////
<<<<<<< HEAD
try:
    run = 1;
    if run:
	#fid_out.write("# Time(ps)\tAverage lipid per area\n");
	psf = '{0}{1}'.format(base_path, structure_file);
    
	flg_top = 0;	    # check bi/mono layer
	flg_btm = 0;	    # check bi/mono layer
    
	timeStamp = [];         # time stamp for trajectory
    
	# data based on trajectory
	STMP  = [];
	#creating BIN 0~90 degree
	BIN = [];
	DIST = [];
	topDIST = [];
	btmDIST = [];
	for ibin in stanalyzer.frange(0, 90, 0.1):
	    BIN.append(ibin);
	    DIST.append(0.0);
	    topDIST.append(0.0);
	    btmDIST.append(0.0);
    
	trj_cnt = 0;
	for idx in range(len(trajectoryFile)):
	    print "[{0}/{1}] trajectory is processing...".format(idx+1, len(trajectoryFile));
	    cnt = 0;
	    # turning on periodic boundary conditions
	    MDAnalysis.core.flags['use_periodic_selections'] = True
	    MDAnalysis.core.flags['use_KDTree_routines'] = False
	    
	    # reading trajectory
	    trj = '{0}{1}'.format(base_path, trajectoryFile[idx]);
	    print 'Reading PSF: ' + psf
	    print 'Reading DCD: ' + trj
	    u = Universe(psf, trj);
	    
	    # read based on frame
	    for ts in u.trajectory:
		# turning on periodic boundary conditions
		MEMB = u.selectAtoms(selQry);
		u.atoms.translate(-MEMB.centerOfMass());
		cnt = cnt + 1;
		trj_cnt = trj_cnt + 1;
		#print "[{0}/{1}] processing...".format(cnt, len(u.trajectory));
		if (cnt % frmInt) == 0:
		    tmp_time = float(trj_cnt) * float(num_ps) - float(num_ps);
		    STMP.append(tmp_time);
		    CRDs = MEMB.coordinates();
		    # for bilayer
		    if len(CRDs) > 0:
			AtRes = getCRDsWithResid(MEMB);
			ResIDs  = AtRes[0]; # residue IDs [resid1, resid2, ..., residN]
			AtmName = AtRes[1]; # [[name11, name12], [name21, name22],....,[nameN1, nameN2]] - assume each residue has 2 atoms
			AtmCrds = AtRes[2] # [[[x11,y11,z11][x12,y12,z12]],[[x21,y21,z21][x22,y22,z22]] ...[[xN1,yN1,zN1][xN2,yN2,zN2]]] 
			for ridx in range(len(ResIDs)):
			    vtx1 = AtmCrds[ridx][0];	# atom C3
			    vtx2 = AtmCrds[ridx][1];	# atom C17
			    #v0 = [0, 0, 1];		# bilayer normal
			    v1 = [];
			    for axi in range(len(vtx1)):
				tmp = vtx1[axi] - vtx2[axi];
				v1.append(tmp);
				
			    if v1[2] < 0:
				v0 = [0, 0, -1];		# bilayer normal
			    else:
				v0 = [0, 0, 1];
    
			    T = toDegree(math.acos(getCosT(v0, v1)));	# degree
			    pos = bisect_left(BIN, T);
			    DIST[pos] = DIST[pos] + 1;
				
		    #print "+----------------------------------------------------------+";
	    
	
	# Writing final output
	all_dist_file = '{0}/all_dist_{1}'.format(out_dir, outFile);
	fid_all_dist = open(all_dist_file, 'w');
	
	# pick one of output file names for DB usage
	out_file = fid_all_dist;
	
	# for all Distribution
	# normalizing DIST
	tmp = [];
	sum_dist = sum(DIST);
	for val in DIST:
	    tmp.append(val/sum_dist);
	DIST = tmp;
	
	cmt_Str = "# Degree\tNumber_of_Sterols\n";
	fid_all_dist.write(cmt_Str);
	for ibin in range(len(BIN)):
	    allStr= "{0}\t{1}\n".format(BIN[ibin], DIST[ibin]);
	    fid_all_dist.write(allStr);
	fid_all_dist.close();
	
	# -------- Drawing graphs
	# Writing Gnuplot script
	outScr = '{0}/gplot{1}.p'.format(out_dir, para_idx);
	outImg  = '{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
	imgPath = "{0}/{1}".format(out_dir, outImg);
	fid_out = open(outScr, 'w');
	gScript = "set terminal png\n";
	gScript = gScript + "set encoding iso_8859_1\n";
	gScript = gScript + "set xlabel 'Degree'\n";
	gScript = gScript + "set ylabel 'Probability'\n";
	gScript = gScript + "set title 'Sterol Tilt'\n";
	gScript = gScript + "set output '{0}'\n".format(imgPath);
	gScript = gScript + """plot "{0}" using 1:2 title "Sterol Tail Tilt" with lines lw 3\n""".format(all_dist_file);
	fid_out.write(gScript);
	fid_out.close();
	
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

=======
run = 1;
if run:
    #fid_out.write("# Time(ps)\tAverage lipid per area\n");
    psf = '{0}{1}'.format(base_path, structure_file);

    flg_top = 0;	    # check bi/mono layer
    flg_btm = 0;	    # check bi/mono layer

    timeStamp = [];         # time stamp for trajectory

    # data based on trajectory
    STMP  = [];
    #creating BIN 0~90 degree
    BIN = [];
    DIST = [];
    topDIST = [];
    btmDIST = [];
    for ibin in stanalyzer.frange(0, 90, 0.1):
	BIN.append(ibin);
	DIST.append(0.0);
	topDIST.append(0.0);
	btmDIST.append(0.0);

    trj_cnt = 0;
    for idx in range(len(trajectoryFile)):
	print "[{0}/{1}] trajectory is processing...".format(idx+1, len(trajectoryFile));
	cnt = 0;
	# turning on periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
	
	# reading trajectory
	trj = '{0}{1}'.format(base_path, trajectoryFile[idx]);
	print 'Reading PSF: ' + psf
	print 'Reading DCD: ' + trj
	u = Universe(psf, trj);
	
	# read based on frame
	for ts in u.trajectory:
	    # turning on periodic boundary conditions
	    MEMB = u.selectAtoms(selQry);
	    u.atoms.translate(-MEMB.centerOfMass());
            cnt = cnt + 1;
	    trj_cnt = trj_cnt + 1;
	    #print "[{0}/{1}] processing...".format(cnt, len(u.trajectory));
	    if (cnt % frmInt) == 0:
		tmp_time = float(trj_cnt) * float(num_ps) - float(num_ps);
		STMP.append(tmp_time);
		CRDs = MEMB.coordinates();
		# for bilayer
		if len(CRDs) > 0:
		    AtRes = getCRDsWithResid(MEMB);
		    ResIDs  = AtRes[0]; # residue IDs [resid1, resid2, ..., residN]
		    AtmName = AtRes[1]; # [[name11, name12], [name21, name22],....,[nameN1, nameN2]] - assume each residue has 2 atoms
		    AtmCrds = AtRes[2] # [[[x11,y11,z11][x12,y12,z12]],[[x21,y21,z21][x22,y22,z22]] ...[[xN1,yN1,zN1][xN2,yN2,zN2]]] 
		    for ridx in range(len(ResIDs)):
			vtx1 = AtmCrds[ridx][0];	# atom C3
			vtx2 = AtmCrds[ridx][1];	# atom C17
			#v0 = [0, 0, 1];		# bilayer normal
			v1 = [];
			for axi in range(len(vtx1)):
			    tmp = vtx1[axi] - vtx2[axi];
			    v1.append(tmp);
			    
			if v1[2] < 0:
			    v0 = [0, 0, -1];		# bilayer normal
			else:
			    v0 = [0, 0, 1];

			T = toDegree(math.acos(getCosT(v0, v1)));	# degree
			pos = bisect_left(BIN, T);
			DIST[pos] = DIST[pos] + 1;
			    
		#print "+----------------------------------------------------------+";
	
    
    # Writing final output
    all_dist_file = '{0}/all_dist_{1}'.format(out_dir, outFile);
    fid_all_dist = open(all_dist_file, 'w');
    
    # pick one of output file names for DB usage
    out_file = fid_all_dist;
    
    # for all Distribution
    # normalizing DIST
    tmp = [];
    sum_dist = sum(DIST);
    for val in DIST:
	tmp.append(val/sum_dist);
    DIST = tmp;
    
    cmt_Str = "# Degree\tNumber_of_Sterols\n";
    fid_all_dist.write(cmt_Str);
    for ibin in range(len(BIN)):
	allStr= "{0}\t{1}\n".format(BIN[ibin], DIST[ibin]);
	fid_all_dist.write(allStr);
    fid_all_dist.close();
    
    # -------- Drawing graphs
    # Writing Gnuplot script
    outScr = '{0}/gplot{1}.p'.format(out_dir, para_idx);
    outImg  = '{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
    imgPath = "{0}/{1}".format(out_dir, outImg);
    fid_out = open(outScr, 'w');
    gScript = "set terminal png\n";
    gScript = gScript + "set encoding iso_8859_1\n";
    gScript = gScript + "set xlabel 'Degree'\n";
    gScript = gScript + "set ylabel 'Probability'\n";
    gScript = gScript + "set title 'Sterol Tilt'\n";
    gScript = gScript + "set output '{0}'\n".format(imgPath);
    gScript = gScript + """plot "{0}" using 1:2 title "Sterol Tail Tilt" with lines lw 3\n""".format(all_dist_file);
    fid_out.write(gScript);
    fid_out.close();
    
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

try:
    print "okay!";
>>>>>>> aa05be30ce412a3a250b73cced1ef91bb83eed20
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

