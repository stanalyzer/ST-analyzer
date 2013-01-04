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

#********************************************
# Define global variables
#********************************************
dbName = 'stanalyzer.db';

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
    def __init__(self, psf, trj):
	self.u = Universe(psf, trj);
	# total number of frames at each trajectory
	self.num_frm = len(self.u.trajectory);
	# total number of atoms
	self.num_atom  = len(self.u.atoms);
	# list of segments
	self.CsegList = self.u.segments;
	self.segList = [];
	for i in range(len(self.CsegList)):
	    self.segList.append(self.CsegList[i].name);
	
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
    funcName = funcName[0:len(funcName)-3];
    fName = [];
    pInfo = [];
    paras = [];
    para_pkey = [];
    print dic["funcName"]
    print funcName
    if funcName in dic["funcName"]:
	print "Function: "
	idx = dic["funcName"].index(funcName);
	print "*** PARA_PKEYS: "
	print dic["para_pkeys"]
	print "({}, {})".format(idx, funcName);
	print "PARMETERS"
	for j in range(len(dic["paraInfo"][idx])):
	    pInfo.append(dic["paraInfo"][idx][j]);
	    paras.append(dic["Paras"][idx][j]);
	    para_pkey.append(dic["para_pkeys"][idx][j]);
	    print " - {}={}".format(dic["paraInfo"][idx][j], dic["Paras"][idx][j]);
	fName = [idx, funcName];
    
    return ([fName, pInfo, paras, para_pkey]);
	