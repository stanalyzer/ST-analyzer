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
num_atoms = paras[4][para_idx];			# pInfo[4] : Total number of atoms
num_atoms = int(num_atoms);
segID     = paras[5][para_idx];
resName   = paras[6][para_idx];
charmm_query = paras[7][para_idx];
cntQry = paras[8][para_idx];			# pInfo[8] : Centering Query
cntAxs = paras[9][para_idx];			# pInfo[9]: Centering axis
#print cntQry;
#print cntAxs

# Reading a trajectory to get the information about atoms
trj = '{0}{1}'.format(base_path, trajectoryFile[0]);
psf = '{0}{1}'.format(base_path, structure_file);

#print "Reading trajectory..."
u = Universe(psf, trj);

#print "creating query..."
selQry = "segid {0} and resname {1} and {2}".format(segID, resName, charmm_query);

#print "Query: {0}".format(selQry);
selAtoms = u.selectAtoms(selQry);

Names = selAtoms.names();					# list names of atoms based on the query
setNames = set(Names);						# make unique list of atoms
listNames = list(setNames);					# convert set to list
sortlistNames = stanalyzer.sort_str_num(listNames, "asc")	# sort them ascending order

num_selatoms = len(Names);
num_uq_selatoms = len(sortlistNames);
num_tails = num_selatoms / num_uq_selatoms;

#print "total {0} tails: {1} atoms and {2} unique atoms".format(num_tails, num_selatoms, num_uq_selatoms);
#print sortlistNames


# defining the order of lipid tail (e.g. CA1, CB1, CC1, etc. Vs. C1, C2)
tmp_num = filter(str.isdigit, sortlistNames[0][1]);
if len(tmp_num) > 0:
    dnst_min = int(sortlistNames[0][1:2]); # reading the first sroted atom and extract the first digit
else:
    dnst_min = filter(str.isdigit, sortlistNames[0]); # reading the first sroted atom and extract the first digit
#print dnst_min

tmp_num = filter(str.isdigit, sortlistNames[len(sortlistNames)-1][1]);
if len(tmp_num) > 0:
    dnst_max = int(sortlistNames[len(sortlistNames)-1][2:]);
else:
    dnst_max = filter(str.isdigit, sortlistNames[len(sortlistNames)-1]);
#print dnst_max

dnst_bin = 1;

#print "dnst_min={}, dnst_max={}, dnst_bin={}".format(dnst_min, dnst_max, dnst_bin);

# query selecting carbon tail
#selQryC2 = "resname DPPE and name C2* and not name C21"; #SCD2
#selQryC3 = "resname DXPE and name C3* and not name C31"; #SCD1
#newQryC2 = "resname DPPE and (name C2* or name H*) and not name C21"; #SCD2
#newQryC3 = "resname DXPE and (name C3* or name H*) and not name C31"; #SCD1

print "ordpara_charmm.py is running...";
#print "AXIS = {}, {}".format(taxis, type(taxis));
#print "Total # atoms = {}, {}".format(num_atoms, type(num_atoms));

#///////////////////////////////////////////////////////////////////////////
# Running actual job
#///////////////////////////////////////////////////////////////////////////

try:
    run = 1;
    if run:
	#out_file = '{0}/{1}'.format(out_dir, outFile);
	#fid_out = open(out_file, 'w')
    
	# describing BIN
	cmt1 = "# Bin ranges\n#"
	#fid_out.write(cmt);
	
	# calculating BIN range
	BIN = [];
	cmt2= '';
	for ibin in frange(dnst_min, dnst_max, dnst_bin):
	    BIN.append(ibin);
	    cmt2= "{0}{1}\t".format(cmt2, ibin);
	    #fid_out.write(cmt);
	cmt2= "{0}\n".format(cmt2);
	
	cmt = "{0}{1}# Range\tDensity\n".format(cmt1, cmt2);
	
	psf = '{0}{1}'.format(base_path, structure_file);
    
	cnt = 0;
	timeStamp = [];         # time stamp for trajectory
	
	# data based on trajectory
	DNST = [];
	STMP = [];
	TOPL = [];	    # order parameter for top layer
	BTML = [];	    # order parameter for bottom layer
	for ibin in frange(dnst_min, dnst_max, dnst_bin):
	    DNST.append(0.0);
	    TOPL.append(0.0);
	    BTML.append(0.0);
    
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
	    allAtoms = u.selectAtoms('all');	# gathering all atoms
	    allCRDs  = allAtoms.coordinates();
	    print '{0} is done!'.format(idx);
	    
	    # read based on frame
	    for ts in u.trajectory:
		#======= Centeralization =========
		if (cntQry != 'no') :
		    #print "Centeralization..."
		    #stanalyzer.centerByCOM(ts, u, cntQry);
		    stanalyzer.centerByRes(ts, u, cntQry, 1, cntAxs); # 1st residue is always chosen for centering membrane
		    #print "DONE!"
		#==================================
		# calculating leaflet
		#L = MDAnalysis.analysis.leaflet.LeafletFinder(u, selQry, cutoff=17.0, pbc=True);
		top_selQry = "{} and (prop z > 0.0)".format(selQry);
		btm_selQry = "{} and (prop z < 0.0)".format(selQry);
		
		#print top_selQry;
		#print btm_selQry;
		topAtoms = u.selectAtoms(top_selQry);
		btmAtoms = u.selectAtoms(btm_selQry);
		
		DegT = [];
		#tclock = cnt;
		cnt = cnt + 1;
		flg_top = 0;	    # set bilayer =1 mono layer = 0
		flg_btm = 0;	    # set bilayer =1 mono layer = 0

		if (cnt % frmInt) == 0:
		    tmp_time = float(cnt) * float(num_ps) - float(num_ps);
		    STMP.append(tmp_time);
		    #print "[{0}ps]selecting atoms...".format(tmp_time);
		    #selAtoms = u.selectAtoms(selQry);
		    #print "DONE!"
		    #if len(L.groups()) > 1:
		    #	topAtoms = L.group(0);		# select atoms in the top layer
		    #--- calculating number ripids
		    if len(topAtoms.resnames()) > 0:
			flg_top = 1;
			
			Names = topAtoms.names();					# list names of atoms based on the query
			setNames = set(Names);						# make unique list of atoms
			listNames = list(setNames);					# convert set to list
			sortlistNames = stanalyzer.sort_str_num(listNames, "asc")	# sort them ascending order
			num_selatoms = len(Names);
			num_uq_selatoms = len(sortlistNames);
			num_top = num_selatoms / num_uq_selatoms;
			#num_top2 = len(topAtoms.resnames()); # length of tail
			#print 'There are {} lipids in top layer'.format(num_top);
			#print 'There are {} lipids in top layer'.format(num_top2);
			
			# find hydrogen atoms based on the current selection
			# get index of current atom selection
			cIndex  = topAtoms.indices();
			for cidx in range(len(cIndex)):
			    h_cnt   = 0;			# hydrogen count
			    tmp_Scd = 0;
			    c_idx = cIndex[cidx]; # current index
			    
			    Cx = allCRDs[c_idx][0];
			    Cy = allCRDs[c_idx][1];
			    Cz = allCRDs[c_idx][2];
			    
			    next_idx = c_idx + 1;
			    
			    while (allAtoms[next_idx].name.upper()[0:1] == 'H'):
				    h_cnt += 1;
				    Hx = allCRDs[next_idx][0];
				    Hy = allCRDs[next_idx][1];
				    Hz = allCRDs[next_idx][2];
				    
				    x = Hx - Cx;
				    y = Hy - Cy;
				    z = Hz - Cz;
				    v1 = [x, y, z];
				    #print "v1:{}, v2:{}".format(v1, v2);
				    cosT = getCosT(v1, v2);
				    tmp_Scd  += (3.0 * cosT * cosT - 1.0);
				    
				    next_idx +=1;
    
			    total_Scd = 0.5 * (tmp_Scd / h_cnt);
			    pos = (cidx) % len(BIN);
			    #print "position={}, h_cnt={}".format(pos, h_cnt);
			    TOPL[pos] += total_Scd;
			    #print TOPL
			    #print "+----------------------------------------------------------+";
		    #btmAtoms = L.group(1);		# select atoms in the bottom layer 
		    if len(btmAtoms.resnames()) > 0:
			flg_btm = 1;
			
			Names = btmAtoms.names();					# list names of atoms based on the query
			setNames = set(Names);						# make unique list of atoms
			listNames = list(setNames);					# convert set to list
			sortlistNames = stanalyzer.sort_str_num(listNames, "asc")	# sort them ascending order
			num_selatoms = len(Names);
			num_uq_selatoms = len(sortlistNames);
			num_btm = num_selatoms / num_uq_selatoms;
			#num_btm2 = len(btmAtoms.resnames()); # length of tail
			#print 'There are {} lipids in bottom layer'.format(num_btm);
			#print 'There are {} lipids in bottom layer'.format(num_btm2);

			# find hydrogen atoms based on the current selection
			# get index of current atom selection
			cIndex  = btmAtoms.indices();
			for cidx in range(len(cIndex)):
			    h_cnt   = 0;			# hydrogen count
			    tmp_Scd = 0;
			    c_idx = cIndex[cidx]; # current index
			    
			    Cx = allCRDs[c_idx][0];
			    Cy = allCRDs[c_idx][1];
			    Cz = allCRDs[c_idx][2];
			    
			    next_idx = c_idx + 1;
			    
			    while (allAtoms[next_idx].name.upper()[0:1] == 'H'):
				    h_cnt += 1;
				    Hx = allCRDs[next_idx][0];
				    Hy = allCRDs[next_idx][1];
				    Hz = allCRDs[next_idx][2];
				    
				    x = Hx - Cx;
				    y = Hy - Cy;
				    z = Hz - Cz;
				    v1 = [x, y, z];
				    #print "v1:{}, v2:{}".format(v1, v2);
				    cosT = getCosT(v1, v2);
				    tmp_Scd  += (3.0 * cosT * cosT - 1.0);
				    
				    next_idx +=1;
    
			    total_Scd = 0.5 * (tmp_Scd / h_cnt);
			    pos = (cidx) % len(BIN);
			    #print "position={}, h_cnt={}".format(pos, h_cnt);
			    BTML[pos] += total_Scd;
			    #print TOPL
			    #print "+----------------------------------------------------------+";
	
    
    
	# Writing final output
	#-----------------------------------------
	#------- This is for top or mono layer
	#-----------------------------------------
	if flg_top > 0:
	    finalTOPL = [];
	    for i in TOPL:
		tmp = abs(i/(len(STMP) * num_top));	# average through trajectory
		finalTOPL.append(tmp);
		
	    # write output into a file
	    top_file = '{0}/top_{1}'.format(out_dir, outFile);
	    fid_top = open(top_file, 'w');
	    for i in range(len(finalTOPL)):
		outStr = "{0}\t{1}\n".format(BIN[i], finalTOPL[i]);
		fid_top.write(outStr);
	    fid_top.close();
	    
	    # -------- Drawing graphs
	    # Writing Gnuplot script
	    outScr = '{0}/top_gplot{1}.p'.format(out_dir, para_idx);
	    outImg  = 'top_{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
	    imgPath = "{0}/{1}".format(out_dir, outImg);
	    fid_out = open(outScr, 'w');
	    gScript = "set terminal png\n";
	    gScript = gScript + "set xlabel 'Carbon Index'\n";
	    gScript = gScript + "set ylabel 'S_CD'\n";
	    gScript = gScript + "set output '{0}'\n".format(imgPath);
	    gScript = gScript + """plot "{0}" using 1:2 title "{1}" with lines lw 3\n""".format(top_file, resName);
	    fid_out.write(gScript);
	    fid_out.close();
	    
	    # Drawing graph with gnuplot
	    subprocess.call(["gnuplot", outScr]);

	    
	#-----------------------------------------
	#------- This is for bottom layer
	#-----------------------------------------
	if flg_btm > 0:
	    finalBTML = [];
	    for i in BTML:
		tmp = abs(i/(len(STMP) * num_btm));	# average through trajectory
		finalBTML.append(tmp);
		
	    # write output into a file
	    btm_file = '{0}/btm_{1}'.format(out_dir, outFile);
	    fid_btm = open(btm_file, 'w');
	    for i in range(len(finalBTML)):
		outStr = "{0}\t{1}\n".format(BIN[i], finalBTML[i]);
		fid_btm.write(outStr);
	    fid_btm.close();
	    
	    # -------- Drawing graphs
	    # Writing Gnuplot script
	    outScr = '{0}/btm_gplot{1}.p'.format(out_dir, para_idx);
	    outImg  = 'btm_{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
	    imgPath = "{0}/{1}".format(out_dir, outImg);
	    fid_out = open(outScr, 'w');
	    gScript = "set terminal png\n";
	    gScript = gScript + "set xlabel 'Carbon Index'\n";
	    gScript = gScript + "set ylabel 'S_CD'\n";
	    gScript = gScript + "set output '{0}'\n".format(imgPath);
	    gScript = gScript + """plot "{0}" using 1:2 title "{1}" with lines lw 3\n""".format(top_file, resName);
	    fid_out.write(gScript);
	    fid_out.close();
	    
	    # Drawing graph with gnuplot
	    subprocess.call(["gnuplot", outScr]);

	#-----------------------------------------
	#------- This is for bilayer
	#-----------------------------------------
	if (flg_top + flg_btm) > 1:
	    # calculating results
	    finalDNST = (np.array(finalTOPL) + np.array(finalBTML)) * 0.5;
	    
	    # write output into a file
	    ave_file = '{0}/ave_{1}'.format(out_dir, outFile);
	    fid_ave = open(ave_file, 'w');
	    for i in range(len(finalDNST)):
		outStr = "{0}\t{1}\n".format(BIN[i], finalDNST[i]);
		fid_ave.write(outStr);
	    fid_ave.close()
	    
	    # Writing Gnuplot script
	    outScr = '{0}/gplot{1}.gpl'.format(out_dir, para_idx);
	    outImg  = 'ave_{0}{1}.png'.format(exe_file[:len(exe_file)-3], para_idx);
	    imgPath = "{0}/{1}".format(out_dir, outImg);
	    fid_out = open(outScr, 'w');
	    gScript = """set terminal png enhanced \n""";
	    gScript = gScript + "set encoding iso_8859_1\n";
	    gScript = gScript + "set output '{0}'\n".format(imgPath);
	    gScript = gScript + "set multiplot layout 3, 1 title 'Order parameter'\n";
	    gScript = gScript + "set tmargin 2\n";
	    
	    gScript = gScript + "set title 'Top Membrane'\n";
	    gScript = gScript + """set xlabel 'Carbon Index'\n""";
	    gScript = gScript + """set ylabel 'S_CD'\n""";
	    gScript = gScript + """plot "{0}" using 1:2 title "{1}" with lines lw 3\n""".format(top_file, resName);
	    
	    gScript = gScript + "set title 'Bottom Membrane'\n";
	    gScript = gScript + "set xlabel 'Carbon Index'\n";
	    gScript = gScript + "set ylabel 'S_CD'\n";
	    gScript = gScript + """plot "{0}" using 1:2 title "{1}" with lines lw 3\n""".format(btm_file, resName);

	    gScript = gScript + "set title 'Both Membrane'\n";
	    gScript = gScript + "set xlabel 'Carbon Index'\n";
	    gScript = gScript + "set ylabel 'S_CD'\n";
	    gScript = gScript + """plot "{0}" using 1:2 title "{1}" with lines lw 3\n""".format(ave_file, resName);
	    
	    gScript = gScript + "unset multiplot\n";
	    gScript = gScript + "set output\n";
	    fid_out.write(gScript);
	    fid_out.close()
	    
	    # Drawing graph with gnuplot
	    subprocess.call(["gnuplot", outScr]);

	# gzip all reaults
	print out_dir;
	print "=== out zip ==="
	outZip = "{0}project_{1}_{2}{3}.tar.gz".format(OUTPUT_HOME, prj_pkey, fName[1], para_idx);
	print outZip
	subprocess.call(["tar", "czf", outZip, out_dir]);
    
	# pick one of output file names for DB usage
	out_file = top_file;
	
	# Update values into gui_outputs
	print "=== DB update ==="
	conn = sqlite3.connect(DB_FILE);
	c    = conn.cursor();
	query = """UPDATE gui_outputs SET status = "Complete", img="{0}", txt="{1}", gzip="{2}" WHERE id = {3}""".format(imgPath, out_file, outZip, pk_output);
	print query
	c.execute(query);
	conn.commit();
	conn.close();
	print "DB update success"
    
    
    
    
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

