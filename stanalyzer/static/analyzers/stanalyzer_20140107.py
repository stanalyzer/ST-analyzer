# for web-interfacer
#import django
#from django.shortcuts		import render_to_response, get_object_or_404, render
#from django.template		import Context, loader, RequestContext
#from django.http		import HttpResponseRedirect, HttpResponse, Http404
from django.utils.encoding	import smart_str, smart_unicode
#from django.utils               import simplejson
#from django.core.urlresolvers	import reverse

# for MDAnalysis
import sys
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import math
import numpy as np

#from MDAnalysis.tests.datafiles import PSF, DCD, PDB, XTC
#import thread
#import Queue
#from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small, PDB, XTC
#import time
#import threading
#from datetime import datetime
# dJango test
#from django import forms
#from django.core.mail import send_mail


# for server side jobs
import os
import string
import random
from random import randint

# for using from object
# https://docs.djangoproject.com/en/dev/topics/forms/?from=olddocs#filefield
#from django import forms
#from django.core.mail import send_mail
#from django.contrib.admin.widgets import FilteredSelectMultiple  
#from django.forms import ModelForm

# importing models
import hashlib
#from django.db import connection, transaction
#from gui.models import User, Project, Job, Parameter

# --- data modification
#cursor.execute("INSERT INTO table_name (column1, column2, column3) VALUES (value1, value2, value3)")
#cursor.execute("DELETE FROM table_name WHERE column1='xxx'")
#cursor.execute("UPDATE table_name SET column1='update it' WHERE column1='xxx' AND column='yyy'")
#transaction.commit_unless_managed()
# --- Data retrieval
#cursor.execute("SELECT column1, column2 FROM table_name WHERE column3 = %s", zzz)
#row = cursor.fetchone()

#for using Ajax: Json
import json
from django.utils import simplejson
import pickle

# for using sqlite3
import sqlite3

# import others
import re
from datetime import datetime
from collections import defaultdict
from bisect import bisect_left

#********************************************
# Define global variables
#********************************************
# name of database
dbName = 'stanalyzer.db';

# dictionary of atomic number
ATOMNUM1 = {'H':1,  'C':6, 'N':7, 'O':8, 'P':15,
	    'S':16, 'K':19 };

ATOMNUM2 = {'LI':3,  'NA':11, 'MG':12, 'CL':17, 'ZN':30,
	    'RB':37, 'CD':48, 'CS':55, 'BA':56, 'HT':1,
	    'HX':1,  'OT':8,  'OX':8, };

ATOMNUM3 = {'LIT':3,  'SOD':11, 'POT':15, 'CLA':20, 'RUB':37,
	    'CAD':48, 'CES':55, 'BAR':56, };

#*******************************************************************
# mapping atoms to thier atomic number in periodic table elements
#*******************************************************************
def firstNletter(N, myStr):
    return myStr[0:N];

def getNletters(N, myList):
    newList = [firstNletter(N, myStr) for myStr in myList];
    return newList;
    
def getAtomNumber(atomList):
    # final atom values corresponding to atomList;
    tmp_atomValues = [0] * len(atomList);
    atomValues = np.array(tmp_atomValues);
    
    # calculating unique atoms     
    uqAtoms = list(set(atomList));
    uqAtomNum = [0] * len(uqAtoms);	# initializing values for uqAtoms
    
    # get atom names with length 3
    uqAtomsL3 = getNletters(3, uqAtoms);
    
    # get atom names with length 2
    uqAtomsL2 = getNletters(2, uqAtoms);

    # get atom names with length 2
    uqAtomsL1 = getNletters(1, uqAtoms);
    
    # find matching atoms listed in the ATOMNUM dictionary
    for key in ATOMNUM3:
	if (key in uqAtomsL3):
	    Idx = findIndex(uqAtomsL3, key);
	    for i in Idx:
		 if uqAtomNum[i] == 0:
		    uqAtomNum[i] = ATOMNUM3[key];

    for key in ATOMNUM2:
	if (key in uqAtomsL2):
	    Idx = findIndex(uqAtomsL2, key);
	    for i in Idx:
		 if uqAtomNum[i] == 0:
		    uqAtomNum[i] = ATOMNUM2[key];

    for key in ATOMNUM1:
	if (key in uqAtomsL1):
	    Idx = findIndex(uqAtomsL1, key);
	    for i in Idx:
		 if uqAtomNum[i] == 0:
		    uqAtomNum[i] = ATOMNUM1[key];
    
    # mapping all atoms to corresponding atom number
    idx = 0;
    for atoms in uqAtoms:
	Idx = findIndex(atomList, atoms);
	atomValues[Idx] = uqAtomNum[idx];
	idx = idx + 1;
    
    return atomValues;
    
#********************************************
# make n-digit random number 
#********************************************
def rand_N_digits(n):
    range_start = 10**(n-1)
    range_end = (10**n)-1
    return randint(range_start, range_end)

#********************************************
# make random string with n letters
#********************************************
def rand_N_letters(n, chars=string.ascii_uppercase + string.digits + string.ascii_lowercase):
    return ''.join(random.choice(chars) for x in range(n))

#********************************************
# find matched item from the given list
#********************************************
def findIndex(myList, strItem):
    Idx = [i for i, x in enumerate(myList) if x == strItem];
    return Idx

#********************************************
# Serverside jobs including
# 1) read directories
#********************************************
class serverside:
    def __init__(self, DIR_TRJ):
    	self.path_trj = DIR_TRJ;
    # default trajectory file
    def showdir(self, *args):
	assert os.path.isdir(self.path_trj)
	if len(args) > 0:
	    assert os.path.isdir(str(args[0]))
	    return os.listdir(str(args[0]));
	else:
	    return os.listdir(self.path_trj)

def createChoices(listObj):
    # get the list of files
    form_list = [];
    for pos, question in enumerate(listObj):
	tmp = [pos, question]
	form_list.append(tmp)
    return form_list

def eval_path(path):
    if path[0] != '/':
        path = '/{0}'.format(path);
    if (path[len(path)-1] != '/'):
        path = "{0}/".format(path);
    return path
        
# calculating angle between p1 and p3 with apex at 1
# p1, p2, and p3 should be numpy array
class geomatric:
    def angle(self, p1, p2, p3): 
	assert (p1.size == p2.size) and (p1.size == p3.size) and (p2.size == p3.size)
	normp1 = math.sqrt(np.sum(p1.__pow__(2)));
	normp3 = math.sqrt(np.sum(p3.__pow__(2)));
	theta  = math.acos(np.sum(p1 * p3)/(normp1 * normp3)) * 180 / math.pi;
	if theta > 90.:
	    theta = 180. - theta
	elif theta < -90.:
	    theta = -1 * (180 + theta)
	return theta;

class simulation:
    # test inputs
    # psf = '/home2/jcjeong/project/stanalyzer1/stanalyzer/trajectory/step5_assembly.psf';
    # dcd = '/home2/jcjeong/project/stanalyzer1/stanalyzer/trajectory/step6_1.dcd';
    
    def __init__(self, psf, trj):
	print "*** simulation ****"
	self.u = Universe(psf, trj);
	
	# total number of frames at each trajectory
	#print "num_frm: "
	self.num_frm = len(self.u.trajectory);
	# time unit per frame
	#print self.num_frm;
	
	#print "num_ps: "
	self.num_ps = np.float16(self.u.trajectory.dt);		# use float16 to fit CHARMM code
	#print self.num_ps;
	
	# total number of atoms
	#print "num_atom"
	self.num_atom  = len(self.u.atoms);
	#print self.num_atom;
	
	# list of segments
	#print "SegList"
	self.CsegList = self.u.segments;
	self.segList = [];
	for i in range(len(self.CsegList)):
	    self.segList.append(self.CsegList[i].name);
	#print self.segList;
	
    def get_segname(self, seg_id):
	return self.segList[int(seg_id)]
    
    def get_seg_residues(self, seg_name):
	selQry = 'segid {0}'.format(seg_name);
	Seg = self.u.selectAtoms(selQry);
	return Seg.residues
	
    
#********************************************
# DB retrieval
#********************************************
class getDB:
    row = [];
    def __init__(self, dbName, user):
    	self.dbName = dbName;
        self.user   = user;
    # default trajectory file
    def get_pbs(self, *args):
        conn = sqlite3.connect(self.dbName);
        c = conn.cursor();
        #print "command is {}".format(args[0]);
        if (args[0] == 'get_recent_pbs'):
            query = "SELECT pbs FROM gui_project WHERE user_id='{0}' ORDER BY date DESC limit 1".format(self.user);
            c.execute(query);
            row = c.fetchall();
            pbs = [];
            for i in row:
                pbs.append(i);
            conn.close();
            return pbs
    

#********************************************
# Find parameters in current function
# funcName : function name to call this function argv[0]
# dic 	   : data read from pickle
#********************************************
def get_myfunction (funcName, dic):
    #print "[get_myfunction]";
    funcName = funcName[0:len(funcName)-3];
    fName = [];
    pInfo = [];
    paras = [];
    para_pkey = [];
    #print "dic[funcName] = {}".format(dic["funcName"]);
    #print "funcName = {}".format(funcName);
    if funcName in dic["funcName"]:
	#print "Function: "
	idx = dic["funcName"].index(funcName);
	#print "*** PARA_PKEYS: "
	#print dic["para_pkeys"]
	#print "(idx={}, funcName={})".format(idx, funcName);
	#print "PARMETERS"
	for j in range(len(dic["paraInfo"][idx])):
	    pInfo.append(dic["paraInfo"][idx][j]);
	    paras.append(dic["Paras"][idx][j]);
	    #print " - {}={}".format(dic["paraInfo"][idx][j], dic["Paras"][idx][j]);

	#print "Para Keys : ";
	#print dic["para_pkeys"][idx];
	for k in range(len(dic["para_pkeys"][idx])):
	    para_pkey.append(dic["para_pkeys"][idx][k]);
	    
	fName = [idx, funcName];
    #print "[----- end get_myfunction ----]"
    return ([fName, pInfo, paras, para_pkey]);


#***************************************
# Generating range for floating numbers
#***************************************
def frange (b, e, i):
    # b: begining
    # e: end
    # i: interval
    
    # make sure they are all floating numbers
    b = float(b);
    e = float(e);
    i = float(i);

    cfloat = b;
    while cfloat <= e:
	yield cfloat
	cfloat += i

#***************************************
# Get atoms belonging to residue IDs
#***************************************
def getCRDsWithResid(AtomGroup):
    CRDs = [];
    ATOMs = [];
    resIDs = AtomGroup.resids();
    resNames = AtomGroup.resnames();
    topCRDs = AtomGroup.coordinates();
    for myid in resIDs:
	myAtoms = [];
	myCRDs = [];
	mycnt = 0;
	for idx in range(len(AtomGroup)):
	    if (AtomGroup[idx].resid == myid):
		myCRDs.append(topCRDs[mycnt]);
		myAtoms.append(AtomGroup[idx].name);
	    mycnt = mycnt + 1;
	ATOMs.append(myAtoms);
	CRDs.append(myCRDs);
    resInfo = [resIDs, ATOMs, CRDs, resNames];	
    return resInfo

#***************************************
# Calculating Items in bin
#***************************************
#from collections import defaultdict
#from bisect import bisect_left
def count_intervals2 (seq, intv ):
    len_bin = len(intv);
    BIN = [];
    for i in range(len(intv)):
	BIN.append(0.0);

    for i in range(len(intv)):
	#print "i ={}, BIN={}".format(i, BIN);
	if i == 0:
	    idx = np.where(seq < intv[i]);
	    BIN[i] = BIN[i] + len(idx[0]);
	    
	elif (i >= len_bin - 1):
	    idx = np.where(seq > intv[i]);
	    BIN[i] = BIN[i] + len(idx[0]);
	else:
	    idx = np.where((seq >= intv[i]) & (seq < intv[i+1]));
	    BIN[i] = BIN[i] + len(idx[0]);
	    
    return BIN;

def count_intervals2_mass (seq, mass, intv ):
    len_bin = len(intv);
    BIN = [];
    for i in range(len(intv)):
	BIN.append(0.0);

    for i in range(len(intv)):
	#print "i ={}, BIN={}".format(i, BIN);
	if i == 0:
	    idx = np.where(seq < intv[i]);
	    BIN[i] = mass[idx[0]].sum();
	    
	elif (i >= len_bin - 1):
	    idx = np.where(seq > intv[i]);
	    BIN[i] = mass[idx[0]].sum();
	else:
	    idx = np.where((seq >= intv[i]) & (seq < intv[i+1]));
	    BIN[i] = mass[idx[0]].sum();
	    
    return BIN;

def count_intervals (sequence, intervals):
    count = defaultdict(int)
    intervals.sort()
    for item in sequence:
        pos = bisect_left(intervals, item)
        if pos == len(intervals):
            count[None] += 1
	    print "WARNING: a value has been found outside min/max of system range! please increase min/max of system range!"
        else:
            count[intervals[pos]] += 1
    return count

def count_intervals_mass (sequence, mass, intervals):
    count = defaultdict(int)
    intervals.sort();
    cnt = 0;
    for item in sequence:
        pos = bisect_left(intervals, item)
        if pos == len(intervals):
            count[None] += mass[cnt];
	    print "WARNING: a value has been found outside min/max of system range! please increase min/max of system range!"
        else:
            count[intervals[pos]] += mass[cnt];
	cnt = cnt + 1;
    return count

# convert degree to radian
def toRadian (fromDeg):
    rad = float(fromDeg) * math.pi / 180.0;
    return rad;

# convert radian to degree
def toDegree(fromRad):
    deg = float(fromRad) * 180.0 / math.pi;
    return deg;

# get cosin theta
def getCosT(v1, v2):
    r1 = 0;
    r2 = 0;
    r1r2 = 0;
    for i in range(len(v1)):
	r1 += v1[i] * v1[i];
	r2 += v2[i] * v2[i];
	r1r2 += v1[i] * v2[i];
    r1 = math.sqrt(r1);
    r2 = math.sqrt(r2);
    cosT = r1r2 / (r1 * r2);
    return cosT;


#********************************************************
# *  File sort containing both string and numbers
#********************************************************
def getSortedList(myList, rev):
    print "def:getSortedList"
    #i[0]: use key index, x[1]: use file name for comparision
    idx = [i[0] for i in sorted(enumerate(myList), key=lambda x:x[1], reverse=rev)]; 
    #print "** func: idx"
    #print idx
    sList = sorted(myList, reverse=rev);
    #print "** func: sList"
    #print sList
    return [idx, sList];

def extStrings(myList):
    print "def:extStrings"
    # find every number and find maximum number in a certain range
    fixDigit = 7;
    newList = [];
    for strLine in myList:
	subStr = '';
	subInt = '';
	num_cnt = 0;
	for i in strLine:
	    #print "*** strLine: {}".format(strLine);
	    #print "ord({0})={1}".format(i, ord(i));
	    if (ord(i) > 47) and (ord(i) < 58):
		num_cnt = num_cnt + 1;
		subInt = "{0}{1}".format(subInt, i);
		#print "subInt : {}".format(subInt);
	    else:
		subStr = "{0}{1}".format(subStr, i);
	
	if num_cnt > 0:
	    num_zero = fixDigit - num_cnt;
	    newNum = '0' * num_zero;
	    newNum = "{0}{1}".format(newNum, subInt);
	    subStr = "{0}{1}".format(subStr, newNum);
	    subStr = "{0}{1}".format(subStr, i);
	    subInt = '';
	    num_cnt = 0;

	#print "subStr: {}".format(subStr);		
	newList.append(subStr);
    return newList;

def sort_str_num(myList, order):
    print "def:sort_str_num"
    #print '{0}-{1}'.format(order, order.upper())

    if (order.upper() == 'DESC'):
	order = True;
    else:
	order = False;
    
    # Sort original list and get the index information
    sListInfo = getSortedList(myList, order);		# False:ascending, True: descending
    sIdx  = sListInfo[0];
    sList = sListInfo[1];
    
    # Extend file name by extending digits with fixed length
    new_sList = extStrings(sList);
    #print new_sList;

    # Sort extended file and get the index of sorted one
    new_sListInfo = getSortedList(new_sList, order);		# False:ascending, True: descending
    new_sIdx = new_sListInfo[0];
    new_sList = new_sListInfo[1];
    
    # Using the index to order original file names
    afterSort = [];
    for i in new_sIdx:
	afterSort.append(sList[i]);
    
    return afterSort;

#********************************************************
# *  Centeralization
#********************************************************
def packintobox(ts):
    x = ts._pos;
    L = ts.dimensions[:3];
    x -= numpy.floor(x/L)*L;

#loose boundary +1 Ans differences
def packintobox2(ts, t_axis):
    x = ts._pos;
    L = ts.dimensions[:3];
    if t_axis == 'x':
	x[:,0] -= numpy.floor(x[:,0]/L[0])*L[0];
    elif t_axis == 'y':
	x[:,1] -= numpy.floor(x[:,1]/L[1])*L[1];
    else:
	x[:,2] -= numpy.floor(x[:,2]/L[2])*L[2];

	
#tight boundary exact system size
def packintobox3(ts, t_axis):
    x = ts._pos;
    L = ts.dimensions[:3];

    # tuple results
    tp_pos_idx = np.where(x[:,2] > L[2]);
    tp_neg_idx = np.where(x[:,2] < 0);

    if len(tp_pos_idx[0]) > 0:
       pos_idx = tp_pos_idx[0];
       x[pos_idx,2] = x[pos_idx,2] - L[2];

    if len(tp_neg_idx[0]) > 0:
       neg_idx = tp_neg_idx[0];
       x[neg_idx,2] = x[neg_idx,2] + L[2];


    
def centerByCOM(ts, u, cntQry):
    MEMB = u.selectAtoms(cntQry);
    com_MEMB = MEMB.centerOfMass();
    
    box = ts.dimensions;
    c_box = 0.5 * box[:3];
    
    t = c_box - com_MEMB;
    
    u.atoms.translate(t);
    packintobox(ts);

    MEMB1 = u.selectAtoms(cntQry);
    com_MEMB1 = MEMB1.centerOfMass();
    
    t = np.array([0,0,0]) - com_MEMB1;
    u.atoms.translate(t);


def centerByRes(ts, u, cntQry, ridx, t_axis):
    
    arr_idx = ridx - 1;
    # move targeted segments into the box
    MEMB = u.selectAtoms(cntQry);
    # pick one residue and calculate center of mass
    #com_Res = MEMB.residues[arr_idx].centerOfMass();
    #com_Res = MEMB.residues[arr_idx].centerOfGeometry();
    mx = MEMB.residues[arr_idx].coordinates()[:,0].min();
    my = MEMB.residues[arr_idx].coordinates()[:,1].min();
    mz = MEMB.residues[arr_idx].coordinates()[:,2].min();
    com_Res = [mx, my, mz];
    box = ts.dimensions[:3];
    
    if t_axis == 'x':
	#print "X"
	# get box information
	bx = 0.5 * box[0];
	by = box[1];
	bz = box[2];
	# calculating distance between system and box
	x = bx - com_Res[0];
	y = 0.0;
	z = 0.0;
	t = np.array([x, y, z]);

    elif t_axis == 'y':
	#print "Y"
	# get box information
	bx = box[0];
	by = 0.5 * box[1];
	bz = box[2];
	# calculating distance between system and box
	x = 0.0;
	y = by - com_Res[1];
	z = 0.0;
	t = np.array([x, y, z]);

    else:
	#print "Z"
	# get box information
	bx = box[0];
	by = box[1];
	bz = 0.5 * box[2];
	# calculating distance between system and box
	x = 0.0;
	y = 0.0;
	z =  bz - com_Res[2];
	t = np.array([x, y, z]);

   
    # move system to the COM of Box and then apply PBC
    u.atoms.translate(t);
    packintobox2(ts, t_axis);
    
    # move entire system by locating COM of MEMB = the COM of Box and than apply PBC
    MEMB1 = u.selectAtoms(cntQry);
    com_MEMB1 = MEMB1.centerOfGeometry();
    
    if t_axis == 'x':
	t1 = np.array(bx - [com_MEMB1[0], 0, 0]);
    elif t_axis == 'y':
	t1 = np.array([0, by - com_MEMB1[1], 0]);
    else:
	t1 = np.array([0, 0, bz - com_MEMB1[2]]);
    
    u.atoms.translate(t1);
    packintobox2(ts, t_axis);
    
    # move entire system by locating COM of MEMB = 0;
    MEMB2 = u.selectAtoms(cntQry);
    com_MEMB2 = MEMB2.centerOfGeometry();
    
    if t_axis == 'x':
	t2 = np.array([com_MEMB2[0], 0, 0]);
    elif t_axis == 'y':
	t2 = np.array([0, com_MEMB2[1], 0]);
    else:
	t2 = np.array([0, 0, com_MEMB2[2]]);

    u.atoms.translate(-t2);

# for trajectories produced by harmm which are already centered
def zeroCenter(ts, u):
    L = ts.dimensions[:3] * 0.5;
    u.atoms.translate(-L);
