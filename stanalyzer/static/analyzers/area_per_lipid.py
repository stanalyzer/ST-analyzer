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
num_atoms = paras[3][para_idx];			# pInfo[3] : Total number of atoms
num_atoms = int(num_atoms);
selQry    = paras[4][para_idx];
sysX   	  = paras[5][para_idx];
sysX 	  = float(sysX);
#sysX 	  = 0.5 * float(sysX);
sysY   	  = paras[6][para_idx];
sysY 	  = float(sysY);
#sysY 	  = 0.5 * float(sysY);
print sysX
print sysY
qull   = paras[7][para_idx];			# pInfo[7] : Location of qull
cntQry = paras[8][para_idx];			# pInfo[8] : Centering Query
cntAxs = paras[9][para_idx];			# pInfo[9]: Centering axis
print cntQry;
print cntAxs


print "Okay I am in area_per_lipid.py!!!!!";
print "Total # atoms = {}, {}".format(num_atoms, type(num_atoms));

#///////////////////////////////////////////////////////////////////////////
# Running actual job
#///////////////////////////////////////////////////////////////////////////
try:
    run = 1;
    if run:
	#fid_out.write("# Time(ps)\tAverage lipid per area\n");
	psf = '{0}{1}'.format(base_path, structure_file);
    
	timeStamp = [];         # time stamp for trajectory
    
	# data based on trajectory
	flg_top = 0;	    # check bi/mono layer
	flg_btm = 0;	    # check bi/mono layer
	STMP  = [];
	TOP_AREA = [];	    # area of individual lipids at top layer
	BTM_AREA = [];	    # area of individual lipids at bottom layer
	TOP_AREA_AVE = [];	    # average area based on residue name at top layer
	BTM_AREA_AVE = [];	    # average area based on residue name at bottom layer
	trj_cnt = 0;
	for idx in range(len(trajectoryFile)):
	    #print "[{0}/{1}] trajectory is processing...".format(idx+1, len(trajectoryFile));
	    cnt = 0;
	    # turning on periodic boundary conditions
	    MDAnalysis.core.flags['use_periodic_selections'] = True
	    MDAnalysis.core.flags['use_KDTree_routines'] = False
	    
	    # reading trajectory
	    trj = '{0}{1}'.format(base_path, trajectoryFile[idx]);
	    #print 'Reading PSF: ' + psf
	    #print 'Reading DCD: ' + trj
	    u = Universe(psf, trj);
	    
	    # read based on frame
	    for ts in u.trajectory:
		
		# reading system size at each frame
		sysX = ts.dimensions[0];
		sysY = ts.dimensions[1];
		sysZ = ts.dimensions[2];
		
		# turning on periodic boundary conditions
		#======= Centeralization =========
		if (cntQry != 'no') :
		    #print "Centeralization..."
		    #stanalyzer.centerByCOM(ts, u, cntQry);
		    stanalyzer.centerByRes2(ts, u, cntQry, 1, cntAxs); # 1st residue is always chosen for centering membrane
		    #print "DONE!"
		#==================================
		
		top_selQry = "{} and (prop z > 0.0)".format(selQry);
		btm_selQry = "{} and (prop z < 0.0)".format(selQry);
		
		#print top_selQry;
		#print btm_selQry;
		topAtoms = u.selectAtoms(top_selQry);
		btmAtoms = u.selectAtoms(btm_selQry);

		cnt = cnt + 1;
		trj_cnt = trj_cnt + 1;
		print "[{0}/{1}] processing...".format(cnt, len(u.trajectory));
		if (cnt % frmInt) == 0: #and cnt == 25:
		    tmp_time = float(trj_cnt) * float(num_ps) - float(num_ps);
		    STMP.append(tmp_time);
		    
		    if len(topAtoms.resnames()) > 0:
			print "TOP: {} residues found!".format(len(topAtoms.resnames()));
			flg_top = 1;
			topArea 	= lipidArea.voroArea(topAtoms, selQry, sysX, sysY); 	# topArea = [resid1,...,residn];
			top_resIDs  = topAtoms.resids();		               		# obtaining residue ids
			top_resNames= topAtoms.resnames();
			top_uq_resNames = lipidArea.stateLipidArea(topAtoms, topArea)[0];	# unique residue to get Average on each Residue [Res1, res2, ..., resN]
			tmp_aveRes  = lipidArea.stateLipidArea(topAtoms, topArea)[1];		# Average area corresponding to uq_resNames [[ave, std, stderr]]
			TOP_AREA.append(topArea);						# TOP_AREA = [topArea_1,...,topArea_n]
			TOP_AREA_AVE.append(tmp_aveRes);					# TOP_AREA_AVE = [[ave1, std1, stderr1], [ave2, std2, stderr2], ...[aveN, stdN, stderrN]]	
			    
		    if len(btmAtoms.resnames()) > 0:
			print "BOTTOM: {} residues found!".format(len(btmAtoms.resnames()));
			flg_btm = 1;
			btmArea 	= lipidArea.voroArea(btmAtoms, selQry, sysX, sysY);	# btmArea = [resid1,...,residn];
			btm_resIDs  = btmAtoms.resids();               			# obtaining residue ids
			btm_resNames= btmAtoms.resnames();
			btm_uq_resNames = lipidArea.stateLipidArea(btmAtoms, btmArea)[0];	# unique residue to get Average on each Residue [Res1, res2, ..., resN]
			tmp_aveRes  = lipidArea.stateLipidArea(btmAtoms, btmArea)[1];	# Average area corresponding to uq_resNames [[ave, std, stderr]]
			BTM_AREA.append(btmArea);						# BTM_AREA = [btmArea_1,...,btmArea_n]
			BTM_AREA_AVE.append(tmp_aveRes);					# BTM_AREA_AVE = [[ave1, std1, stderr1], [ave2, std2, stderr2], ...[aveN, stdN, stderrN]]
	    
	
	# Writing final output
	top_res_file = '{0}/top_res_{1}'.format(out_dir, outFile);
	top_ave_file = '{0}/top_ave_{1}'.format(out_dir, outFile);
	btm_res_file = '{0}/btm_res_{1}'.format(out_dir, outFile);
	btm_ave_file = '{0}/btm_ave_{1}'.format(out_dir, outFile);
	fid_top_res = open(top_res_file, 'w');
	fid_top_ave = open(top_ave_file, 'w');
	fid_btm_res = open(btm_res_file, 'w');
	fid_btm_ave = open(btm_ave_file, 'w');
	
	# pick one of output file names for DB usage
	out_file = top_ave_file;
	
	if flg_top > 0:
	    cmt_aveStr = "# Time(ps)";
	    for i in top_uq_resNames:
		cmt_aveStr ="{0}\t{1}".format(cmt_aveStr, i);
	    cmt_aveStr = "{0}\n".format(cmt_aveStr);
    
	    cmt_resStr = "# Time(ps)";
	    for i in top_resIDs:
		cmt_resStr ="{0}\t{1}".format(cmt_resStr, i);
	    cmt_resStr = "{0}\tAVERAGE\n".format(cmt_resStr);
    
	    fid_top_ave.write(cmt_aveStr);
	    fid_top_res.write(cmt_resStr);
	    
	    tstamp = 0;
	    for perRes in TOP_AREA_AVE:
		outStr = "{0}".format(STMP[tstamp]);
		tstamp = tstamp + 1;
		for istat in perRes:
		    outStr = "{0}\t{1}".format(outStr, istat[0]);
		outStr= "{0}\n".format(outStr);
		fid_top_ave.write(outStr);
	    fid_top_ave.close()
	    TOP_AVE_VAL = np.array(TOP_AREA_AVE);
	    max_top = TOP_AVE_VAL[:,0][:,0].max();
	    min_top = TOP_AVE_VAL[:,0][:,0].min();
	    
	    tstamp = 0;
	    for resIdx in TOP_AREA:
		outStr = "{0}".format(STMP[tstamp]);
		tstamp = tstamp + 1;
		for j in range(len(resIdx)):
		    outStr = "{0}\t{1}".format(outStr, resIdx[j]);
		outStr= "{0}\n".format(outStr);
		fid_top_res.write(outStr);
	    fid_top_res.close()
	
	if flg_btm > 0:
	    cmt_aveStr = "# Time(ps)";
	    for i in btm_uq_resNames:
		cmt_aveStr ="{0}\t{1}".format(cmt_aveStr, i);
	    cmt_aveStr = "{0}\n".format(cmt_aveStr);
    
	    cmt_resStr = "# Time(ps)";
	    for i in btm_resIDs:
		cmt_resStr ="{0}\t{1}".format(cmt_resStr, i);
	    cmt_resStr = "{0}\tAVERAGE\n".format(cmt_resStr);
	    
	    fid_btm_ave.write(cmt_aveStr);
	    fid_btm_res.write(cmt_resStr);
	    tstamp = 0;
	    for perRes in BTM_AREA_AVE:
		outStr = "{0}".format(STMP[tstamp]);
		tstamp = tstamp + 1;
		for istat in perRes:
		    outStr = "{0}\t{1}".format(outStr, istat[0]);
		outStr= "{0}\n".format(outStr);
		fid_btm_ave.write(outStr);
	    fid_btm_ave.close()
	    BTM_AVE_VAL = np.array(BTM_AREA_AVE);
	    max_btm = BTM_AVE_VAL[:,0][:,0].max();
	    min_btm = BTM_AVE_VAL[:,0][:,0].min();
	    
	    tstamp = 0;
	    for resIdx in BTM_AREA:
		outStr = "{0}".format(STMP[tstamp]);
		tstamp = tstamp + 1;
		for j in range(len(resIdx)):
		    outStr = "{0}\t{1}".format(outStr, resIdx[j]);
		outStr= "{0}\n".format(outStr);
		fid_btm_res.write(outStr);
	    fid_btm_res.close()
    
    
	# -------- Drawing graphs
	if flg_top > 0:
	    outScr = '{0}/gplot{1}.p'.format(out_dir, para_idx);
	    outImg  = '{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
	    imgPath = "{0}/{1}".format(out_dir, outImg);
	    fid_out = open(outScr, 'w');
	    gScript = """set terminal png enhanced \n""";
	    gScript = gScript + "set output '{0}'\n".format(imgPath);
	    gScript = gScript + "set multiplot layout 2, 1 title 'Area per lipid'\n";
	    
	    # top layer
	    num_top_tics = 3.0;
	    intx = (max_top - min_top) / num_top_tics;
	    # if graph is broken please remove margin
	    #gScript = gScript + """set tmargin at screen 0.93; set bmargin at screen 0.54\n""";
	    #gScript = gScript + """set lmargin at screen 0.20; set rmargin at screen 0.85\n""";
	    gScript = gScript + """set xtics offset 0,0.5; unset xlabel\n""";
	    
	    if intx >= 0.0001:
		gScript = gScript + """set ytics {0:10.4f},{1:10.4f},{2:10.4f}; unset ylabel\n""".format(min_top, intx, max_top);

	    gScript = gScript + """set label 1 'Top' at graph 0.01, 0.95 font ',8'\n""";
	    gScript = gScript + """set ylabel 'Area per Lipid [A^2 ]' offset 0,-8\n""";
	    gScript = gScript + """plot "{0}" using ($1*0.001):2 title "{1}" with lines lw 3""".format(top_ave_file, top_uq_resNames[0]);
	    rcnt = 2;
	    for ridx in range(len(top_uq_resNames)-1):
		rcnt = rcnt + 1;
		gScript = """{0}, "{1}" using  ($1*0.001):{2} title "{3}" with lines lw 3""".format(gScript, top_ave_file, rcnt, top_uq_resNames[ridx+1]);
	    gScript = """{0}\n""".format(gScript);
	    
	if flg_btm > 0:
	    num_btm_tics = 3.0;
	    intx = (max_top - min_btm) / num_btm_tics;
	    # if graph is broken please remove margin
	    #gScript = gScript + """set tmargin at screen 0.48; set bmargin at screen 0.08\n""";
	    #gScript = gScript + """set lmargin at screen 0.20; set rmargin at screen 0.85\n""";
	    gScript = gScript + """set xtics offset 0,0.5; unset xlabel\n""";
	    
	    if intx >= 0.0001:
		gScript = gScript + """set ytics {0:10.4f},{1:10.4f},{2:10.4f}; unset ylabel\n""".format(min_btm, intx, max_btm);
	    
	    gScript = gScript + """set label 1 'Bottom' at graph 0.01, 0.95 font ',8'\n""";
	    gScript = gScript + """set xlabel 'Time (ns)' offset 0,1\n""";
	    gScript = gScript + """plot "{0}" using  ($1*0.001):2 title "{1}" with lines lw 3""".format(btm_ave_file, btm_uq_resNames[0]);
	    rcnt = 2;
	    for ridx in range(len(btm_uq_resNames)-1):
		rcnt = rcnt + 1;
		gScript = """{0}, "{1}" using  ($1*0.001):{2} title "{3}" with lines lw 3""".format(gScript, btm_ave_file, rcnt, btm_uq_resNames[ridx+1]);
	    gScript = """{0}\n""".format(gScript);
	    
	    gScript = gScript + "unset multiplot\n";
	    gScript = gScript + "set output\n";
	    fid_out.write(gScript);
	    fid_out.close()
	    
	if (flg_top + flg_btm) > 1:
	    # Writing Gnuplot script
	    outScr = '{0}/gplot{1}.gpl'.format(out_dir, para_idx);
	    outImg  = '{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
	    imgPath = "{0}/{1}".format(out_dir, outImg);
	    fid_out = open(outScr, 'w');
	    gScript = """set terminal png enhanced \n""";
	    gScript = gScript + "set output '{0}'\n".format(imgPath);
	    gScript = gScript + "set multiplot layout 2, 1 title 'Area per lipid'\n";
	    
	    # top layer
	    num_top_tics = 3.0;
	    intx = (max_top - min_top) / num_top_tics;
	    # if graph is broken please remove margin
	    #gScript = gScript + """set tmargin at screen 0.93; set bmargin at screen 0.54\n""";
	    #gScript = gScript + """set lmargin at screen 0.20; set rmargin at screen 0.85\n""";
	    gScript = gScript + """set xtics offset 0,0.5; unset xlabel\n""";
	    
	    if intx >= 0.0001:
		gScript = gScript + """set ytics {0:10.4f},{1:10.4f},{2:10.4f}; unset ylabel\n""".format(min_top, intx, max_top);
	    
	    gScript = gScript + """set label 1 'Top' at graph 0.01, 0.95 font ',8'\n""";
	    gScript = gScript + """set ylabel 'Area per Lipid [A^2 ]' offset 0,-8\n""";
	    gScript = gScript + """plot "{0}" using ($1*0.001):2 title "{1}" with lines lw 3""".format(top_ave_file, top_uq_resNames[0]);
	    rcnt = 2;
	    for ridx in range(len(top_uq_resNames)-1):
		rcnt = rcnt + 1;
		gScript = """{0}, "{1}" using  ($1*0.001):{2} title "{3}" with lines lw 3""".format(gScript, top_ave_file, rcnt, top_uq_resNames[ridx+1]);
	    gScript = """{0}\n""".format(gScript);
	    
	    # bottom layer
	    num_btm_tics = 3.0;
	    intx = (max_btm - min_btm) / num_btm_tics;
	    # if graph is broken please remove margin
	    #gScript = gScript + """set tmargin at screen 0.48; set bmargin at screen 0.08\n""";
	    #gScript = gScript + """set lmargin at screen 0.20; set rmargin at screen 0.85\n""";
	    gScript = gScript + """set xtics offset 0,0.5; unset xlabel\n""";
	    
	    if intx >= 0.0001:
		gScript = gScript + """set ytics {0:10.4f},{1:10.4f},{2:10.4f}; unset ylabel\n""".format(min_btm, intx, max_btm);
	    
	    gScript = gScript + """set label 1 'Bottom' at graph 0.01, 0.95 font ',8'\n""";
	    gScript = gScript + """set xlabel 'Time (ns)' offset 0,1\n""";
	    gScript = gScript + """plot "{0}" using  ($1*0.001):2 title "{1}" with lines lw 3""".format(btm_ave_file, btm_uq_resNames[0]);
	    rcnt = 2;
	    for ridx in range(len(btm_uq_resNames)-1):
		rcnt = rcnt + 1;
		gScript = """{0}, "{1}" using  ($1*0.001):{2} title "{3}" with lines lw 3""".format(gScript, btm_ave_file, rcnt, btm_uq_resNames[ridx+1]);
	    gScript = """{0}\n""".format(gScript);
	    
	    gScript = gScript + "unset multiplot\n";
	    gScript = gScript + "set output\n";
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

