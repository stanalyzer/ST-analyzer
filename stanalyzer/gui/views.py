#/usr/bin/env python

#----------------------------------------------------------------------------------------
# Author: Jong Cheol Jeong (korjcjeong@yahoo.com, people.eecs.ku.edu/~jjeong)
# 	  Bioinformatics center, The University of Kansas
#----------------------------------------------------------------------------------------

# for web-interfacer
from django.shortcuts		import render_to_response, get_object_or_404, render
from django.template		import Context, loader, RequestContext
from django.http		import HttpResponseRedirect, HttpResponse, Http404
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
import stat
import string
import random
import math
import shutil
from random import randint
from os import path

# for using from object
# https://docs.djangoproject.com/en/dev/topics/forms/?from=olddocs#filefield
#from django import forms
#from django.core.mail import send_mail
#from django.contrib.admin.widgets import FilteredSelectMultiple  
#from django.forms import ModelForm

# importing models
import hashlib
from django.db import connection, transaction
from gui.models import User, Project, Job, Parameter
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

# download
import os, mimetypes
from django.core.servers.basehttp import FileWrapper

#********************************************
# Define global variables
#********************************************
#
#dbName = settings.DATABASES["default"]["NAME"];
SESSION_TIME_OUT = 3000;
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__));
dbName = "{}/stanalyzer.db".format(PROJECT_ROOT[0:len(PROJECT_ROOT)-4]);
MEDIA_HOME = "{}/media".format(PROJECT_ROOT[0:len(PROJECT_ROOT)-4]);
ANALYZER_HOME = "{}/static/analyzers".format(PROJECT_ROOT[0:len(PROJECT_ROOT)-4]);
#********************************************
# Parsing wrapped lists
# This function is used for parsing parameters and thier information
#********************************************
def parseWrapList(strList):
    num_list = len(strList);
    newList = [];
    #print num_list;
    for i in range(num_list):
        #print "Iteration[{}]".format(i);
        strElm = strList[i];
        #print strElm
        lstElm = strElm.split(',');      # split based on comma
        #print lstElm
        newList.append(lstElm);
    return newList;
        

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
def rand_N_letters(n, chars=string.letters + string.digits):
    newStr = ''.join(random.choice(chars) for x in range(n));
    return newStr


#********************************************
# Serverside jobs including
# 1) read directories
#********************************************
class serverside:
    def __init__(self, *args):
        if len(args) > 0:
            self.path_trj = str(args[0]);
    # default trajectory file

    def showdir(self, *args):
        try:
            if len(args) > 0:
                #assert os.path.isdir(str(args[0]))
                return os.listdir(str(args[0]));
            else:
                return os.listdir(self.path_trj)
        except:
            return ['n/a']

    def isvalidate(self, *args):
        status = [];
        try:
            tmp = os.path.isdir(self.path_trj);
            if tmp:
                status = ['dir'];
            tmp = os.path.isfile(self.path_trj);
            if tmp:
                status = ['file'];
        except:
            status =[];
        return status;
    
    def permission_write(self):
        st = os.stat(self.path_trj);
        owner = bool(st.st_mode & stat.S_IWUSR);
        group = bool(st.st_mode & stat.S_IWGRP);
        other = bool(st.st_mode & stat.S_IWOTH);
        return [owner, group, other]

    def permission_exec(self):
        st = os.stat(self.path_trj);
        owner = bool(st.st_mode & stat.S_IXUSR);
        group = bool(st.st_mode & stat.S_IXGRP);
        other = bool(st.st_mode & stat.S_IXOTH);
        return [owner, group, other]

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
# Delete files and directories
#********************************************
def delDir(dir_path):
    print "+++ delDir"
    if os.path.isdir(dir_path):
	for f in os.listdir(dir_path):
	    file_path = os.path.join(dir_path, f)
	    try:
		if os.path.isfile(file_path):
		    #os.remove(file_path);
		    print "remove {}".format(file_path);
		else:
		    shutil.rmtree(file_path);
	    except Exception, e:
		print e
	#shutil.rmtree(dir_path);
	print "remove {}".format(dir_path);


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

def getUsrLevel (dbName, user_id):
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    query = "SELECT level FROM gui_user WHERE uid='{0}'".format(user_id);
    c.execute(query);
    row = c.fetchone();
    level = row[0];
    conn.close();
    return level

def delOutputs (IDs, dbName, delFlg):
    print "*** FUNC: delOutputs"
    print "FUNC: delOutputs"
    # delete gui_job
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    if isinstance(IDs, list):
	for output_id in IDs:
	    if delFlg == 'true':
		query = "SELECT img, txt, gzip FROM gui_outputs WHERE id={0}".format(output_id);
		c.execute(query);
		row = c.fetchall();
		for f in row:
		    if '/' in f[0]:
			path_img = f[0].split('/');
			path_dir = '/'.join(path_img[:len(path_img)-1])
			delDir(path_dir);

	    query = "DELETE FROM gui_outputs WHERE id={0}".format(output_id);
	    print query;
    else:
	if delFlg == 'true':
	    query = "SELECT img, txt, gzip FROM gui_outputs WHERE id={0}".format(IDs);
	    c.execute(query);
	    row = c.fetchall();
	    for f in row:
		if '/' in f[0]:
		    path_img = f[0].split('/');
		    path_dir = '/'.join(path_img[:len(path_img)-1])
		    delDir(path_dir);

	query = "DELETE FROM gui_outputs WHERE id={0}".format(IDs);
	print query;
    conn.close();

def delOutputs_from_jobID (IDs, dbName, delFlg):
    print "*** FUNC: delOutputs_from_jobID"
    # delete gui_job
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    if isinstance(IDs, list):
	print "--- this runs with LIST"
	for job_id in IDs:
	    if delFlg == 'true':
		query = "SELECT img, txt, gzip FROM gui_outputs WHERE job_id={0}".format(job_id);
		c.execute(query);
		row = c.fetchall();
		#print "====== Split word ===="
		for f in row:
		    if '/' in f[0]:
			path_img = f[0].split('/');
			path_dir = '/'.join(path_img[:len(path_img)-1])
			delDir(path_dir);

	    query = "DELETE FROM gui_outputs WHERE job_id={0}".format(job_id);
	    print query;
    else:
	print "--- this runs with SCHOLAR {0}".format(IDs);
	print delFlg
	print "Type is {0}".format(type(delFlg));
	if delFlg == 'true':
	    query = "SELECT img, txt, gzip FROM gui_outputs WHERE job_id={0}".format(IDs);
	    print query
	    c.execute(query);
	    row = c.fetchall();
	    print row
	    for f in row:
		print "====== Split word ===="
		if '/' in f[0]:
		    path_img = f[0].split('/');
		    path_dir = '/'.join(path_img[:len(path_img)-1])
		    print path_dir
		    delDir(path_dir);
		else:
		    print "With no path {0}".format(f[0]);
		    
	query = "DELETE FROM gui_outputs WHERE job_id={0}".format(IDs);
	print query;
    conn.close();

def delJobs (IDs, dbName, delFlg):
    print "*** FUNC: delJobs"
    # delete gui_job
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    if isinstance(IDs, list):
	for job_id in IDs:
	    #delete gui_parameter, gui_outputs
	    query = "DELETE FROM gui_parameter WHERE job_id={0}".format(job_id);
	    print query
	    delOutputs_from_jobID(job_id, dbName, delFlg);
	    if delFlg == 'true':
		#deleting job output directory
		query = "SELECT output FROM gui_job WHERE id={0}".format(job_id);
		c.execute(query);
		row = c.fetchall();
		for f in row:
		    if '/' in f[0]:
			delDir(f[0]);
	    query = "DELETE FROM gui_job WHERE id={0}".format(job_id);
	    print query
    else:
	#delete gui_parameter, gui_outputs
	query = "DELETE FROM gui_parameter WHERE job_id={0}".format(IDs);
	print query
	delOutputs_from_jobID(IDs, dbName, delFlg);
	if delFlg == 'true':
	    query = "SELECT output FROM gui_job WHERE id={0}".format(IDs);
	    c.execute(query);
	    row = c.fetchall();
	    for f in row:
		if '/' in f[0]:
		    delDir(f[0]);
	query = "DELETE FROM gui_job WHERE id={0}".format(IDs);
	print query
    
    conn.close();
    
def delJobs_from_prjID (IDs, dbName, delFlg):
    print "*** FUNC: delJobs_from_prjID"
    # delete gui_job
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    if isinstance(IDs, list):
	for proj_id in IDs:
	    query = "SELECT id FROM gui_job WHERE proj_id = {0}".format(proj_id);
	    c.execute(query);
	    row = c.fetchall();
	    for item in row:
		job_id = item[0];
		#delete gui_parameter, gui_outputs
		query = "DELETE FROM gui_parameter WHERE job_id={0}".format(job_id);
		print query
		delOutputs_from_jobID(job_id, dbName, delFlg);
		if delFlg == 'true':
		    query = "SELECT output FROM gui_job WHERE id={0}".format(job_id);
		    print query
		    c.execute(query);
		    row2 = c.fetchall();
		    print "Print Row2"
		    print row2
		    for f in row2:
			if '/' in f[0]:
			    delDir(f[0]);

	    # delte gui_job corresponding to proj_id
	    query = "DELETE FROM gui_job WHERE proj_id={0}".format(proj_id);
	    print query
    else:
	query = "SELECT id FROM gui_job WHERE proj_id = {0}".format(IDs);
	c.execute(query);
	row = c.fetchall();
	for item in row:
	    job_id = item[0];
	    #delete gui_parameter, gui_outputs
	    query = "DELETE FROM gui_parameter WHERE job_id={0}".format(job_id);
	    print query
	    delOutputs_from_jobID(job_id, dbName, delFlg);
	if delFlg == 'true':
	    query = "SELECT output FROM gui_job WHERE id={0}".format(IDs);
	    c.execute(query);
	    row2 = c.fetchall();
	    for f in row2:
		if '/' in f[0]:
		    delDir(f[0]);
	# delte gui_job corresponding to proj_id
	query = "DELETE FROM gui_job WHERE proj_id={0}".format(IDs);
	print query
    conn.close();

def delProjects (IDs, dbName, delFlg):
    print "*** FUNC: delProjects"
    delJobs_from_prjID(IDs, dbName, delFlg);
    if isinstance(IDs, list):
	# delete gui_path_input, gui_path_output, gui_path_python
	for proj_id in IDs:
	    query = "DELETE FROM gui_path_input WHERE proj_id = {0}".format(proj_id);
	    print query
	    query = "DELETE FROM gui_path_output WHERE proj_id = {0}".format(proj_id);
	    print query
	    query = "DELETE FROM gui_path_python WHERE proj_id = {0}".format(proj_id);
	    print query
	    query = "DELETE FROM gui_project WHERE id = {0}".format(proj_id);
	    print query
    else:
	query = "DELETE FROM gui_path_input WHERE proj_id = {0}".format(IDs);
	print query
	query = "DELETE FROM gui_path_output WHERE proj_id = {0}".format(IDs);
	print query
	query = "DELETE FROM gui_path_python WHERE proj_id = {0}".format(IDs);
	print query
	query = "DELETE FROM gui_project WHERE id = {0}".format(IDs);
	print query
		

#********************************************************
# *  Validating path
#********************************************************
def pathValidation(request):
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
        if (cmd == 'path_validation'):
            path  = request.POST.get('path');
            #print path
            server = serverside(path);
            status = server.isvalidate();
            #print status;
            #print "Status = {}. len={}".format(status, len(status));
            #if len(status) > 0:
            #    print status[0];
                
        c = {
                  'type'           : status,
            }
        return HttpResponse(json.dumps(c));

#********************************************************
# *  Validating path
#********************************************************
def info_permission_write(request):
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
        if (cmd == 'permission_write'):
            path  = request.POST.get('path');
            #print path
            server = serverside(path);
            status = server.permission_write();
            #print status;                           # status = [owner, group, other]
            #print "Status = {}. len={}".format(status, len(status));
            if len(status) > 0:
                print status[0];
                
        c = {
                  'write_info'           : status,
            }
        return HttpResponse(json.dumps(c));

def info_permission_exec(request):
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
        if (cmd == 'permission_exec'):
            path  = request.POST.get('path');
            #print path
            server = serverside(path);
            status = server.permission_exec();
            #print status;                           # status = [owner, group, other]
            #print "Status = {}. len={}".format(status, len(status));
            if len(status) > 0:
                print status[0];
                
        c = {
                  'exec_info'           : status,
            }
        return HttpResponse(json.dumps(c));

#********************************************************
# *  GET DB Information
#********************************************************
def getDBinfo(request):
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
        pbs = [];
        if (cmd == 'get_recent_pbs'):
            get_db = getDB(dbName, request.session['user_id']);
            pbs = get_db.get_pbs(cmd);
	if len(pbs) < 1:
	    pbs = "#!/bin/csh\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=500mb\n#PBS -l walltime=72:00:00\n#PBS -l cput=72:00:00\n#PBS -q default\n";
	    
        c = {
                  'pbs'           : pbs,
              }
        request.session.set_expiry(SESSION_TIME_OUT);        # session will be closed in 10 minutes
        return HttpResponse(json.dumps(c));


def help(request):
    template = 'gui/help.html';
    c = {
	    'pgTitle'	    : 'ST-Analyzer!',
	}
    return render_to_response(template, c, context_instance = RequestContext(request) )

def index(request):
    template = 'gui/login.html';
    c = {
	    'pgTitle'	    : 'MBanalyzer test page!',
	}
    return render_to_response(template, c, context_instance = RequestContext(request) )

def desktop(request):
    template = 'desktop/index.html';
    c = {
	    'pgTitle'	    : 'MBanalyzer test page!',
	}
    return render_to_response(template, c, context_instance = RequestContext(request) )

#********************************************************
# *  for logout
#********************************************************
def logout(request):
    # initializing all sessions
    for sesskey in request.session.keys():
        del request.session[sesskey]
    request.session.set_expiry(1); 
    return HttpResponseRedirect("/gui/")

#********************************************************
# *  for login GUI
#********************************************************
def login(request):
    # initializing all sessions
    for sesskey in request.session.keys():
        del request.session[sesskey];
    
    userid   = request.POST.get('userID');
    pwd      = request.POST.get('pwd');
    userid = userid.strip(' \t\n\r');
    pwd    = pwd.strip(' \t\n\r');
    if (userid == '') or (pwd == ''):
        template = 'gui/login.html';
        Msg      = 'user ID and password should be given!'
    else:
        conn = sqlite3.connect(dbName);
        c = conn.cursor();
        c.execute("select * from gui_user")
        row = c.fetchone();
        #print dbName
        #print row
        
        if not row:
            global user
            user = 'admin';
            pwd   = '12345';
            email = 'admin@stanalyzer.org'
            level = '10'
            hpwd = hashlib.md5(pwd);
            hpwd = hpwd.hexdigest();
            query = "INSERT INTO gui_user (uid, pwd, email, level) VALUES ('{0}', '{1}', '{2}', {3})".format(user, hpwd, email, level)
            #print query
            c.execute(query)
            conn.commit();
            Msg = 'admin has been initialized!'
            template = 'gui/login.html';
        else:
            hpwd = hashlib.md5(pwd);
            hpwd = hpwd.hexdigest();
            query = "SELECT uid, pwd from gui_user WHERE uid='{0}' and pwd='{1}'".format(userid, hpwd);
            #print query
            c.execute(query)
            row = c.fetchone();
            if not row:
                template = 'gui/login.html';
                Msg = 'ID and PWD do not match!';
            else:
                #template = 'gui/stanalyzer.html';
                request.session['user_id'] = userid;
                request.session.set_expiry(SESSION_TIME_OUT);        # session will be closed in 10 minutes
                #Msg = "{0}={1} and {2}={3}".format(userid, row[0], pwd, row[1]);
                Msg = userid
                conn.close();
                return HttpResponseRedirect("/gui/stanalyzer/")
                #return HttpResponseRedirect("/gui/desktop/")
        conn.close();
    c = {
	    'errMsg'	    : Msg,
	}
    return render_to_response(template, c, context_instance = RequestContext(request) )


    
    
#********************************************************
# *  for Project View
#********************************************************
def toyView_data(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )


    request.session.set_expiry(SESSION_TIME_OUT);
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    
    #----------------------------------------
    # INITIALIZING gui_user table
    #----------------------------------------
    c.execute("DELETE FROM gui_user");
    conn.commit();
    
    pwd = "12345";
    hpwd = hashlib.md5(pwd);
    hpwd = hpwd.hexdigest();
    query = """INSERT INTO gui_user (uid, pwd, email, level) \
            VALUES ("admin", "{0}", "admin@stanalyzer.org", 10) """.format(hpwd);
    c.execute(query);
    conn.commit();
    query = """INSERT INTO gui_user (uid, pwd, email, level) \
            VALUES ("jcjeong", "{0}", "jcjeong@stanalyzer.org", 10) """.format(hpwd);
    c.execute(query);
    conn.commit();
    query = """INSERT INTO gui_user (uid, pwd, email, level) \
            VALUES ("guest", "{0}", "guest@stanalyzer.org", 0) """.format(hpwd);
    c.execute(query);
    conn.commit();
    print "gui_user initialized!"
    
    #----------------------------------------
    # INITIALIZING gui_project
    #----------------------------------------
    c.execute("DELETE FROM gui_project");
    conn.commit();
    cdate = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
    pbs = "#!/bin/csh\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=500mb\n#PBS -l walltime=72:00:00\n#PBS -l cput=72:00:00\n#PBS -q default\n";
    query = """INSERT INTO gui_project (user_id, name, date, pbs) \
            VALUES ("admin", "project toy example1", "{0}", "{1}") """.format(cdate, pbs);
    c.execute(query);
    conn.commit();
    print "gui_project initialized!"
    
    #----------------------------------------
    # INITIALIZING gui_path_input
    #----------------------------------------
    c.execute("DELETE FROM gui_path_input");
    conn.commit();
    path = "/home2/jcjeong/project/stanalyzer2/stanalyzer/trajectory";
    query = """INSERT INTO gui_path_input (proj_id, path) \
            VALUES (1, "{0}") """.format(path);
    c.execute(query);
    conn.commit();
    print "gui_path_input initialized!"

    #----------------------------------------
    # INITIALIZING gui_path_input
    #----------------------------------------
    c.execute("DELETE FROM gui_path_input");
    conn.commit();
    path = "/home2/jcjeong/project/stanalyzer2/stanalyzer/trajectory";
    query = """INSERT INTO gui_path_input (proj_id, path) \
            VALUES (1, "{0}") """.format(path);
    c.execute(query);
    conn.commit();
    print "gui_path_input initialized!"

    #----------------------------------------
    # INITIALIZING gui_path_output
    #----------------------------------------
    c.execute("DELETE FROM gui_path_output");
    conn.commit();
    path = "/home2/jcjeong/tmp";
    query = """INSERT INTO gui_path_output (proj_id, path) \
            VALUES (1, "{0}") """.format(path);
    c.execute(query);
    conn.commit();
    print "gui_path_output initialized!"
    
    #----------------------------------------
    # INITIALIZING gui_path_python
    #----------------------------------------
    c.execute("DELETE FROM gui_path_python");
    conn.commit();
    path = "/home/sunhwan/local/python/bin/python";
    query = """INSERT INTO gui_path_python (proj_id, path) \
            VALUES (1, "{0}") """.format(path);
    c.execute(query);
    conn.commit();
    print "gui_path_python initialized!"
    conn.close();
    
    return HttpResponseRedirect("/")

    
def toyView_prj(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
    else:
        cmd = 'http';           # connection with URL

    #---------------------------------------------#
    # HANDLING PROJECT TABLE ** START **
    #---------------------------------------------#
    # Initializing variables
    parents = [];
    indent = [];
    parent = [];
    title  = [];
    path1  = [];
    path2  = [];
    path3  = [];
    path4  = [];
    path5  = [];
    pbs    = [];
    date   = [];
    numRec = [];
    
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    c.execute("select id, user_id, name, path1, path2, path3, path4, path5, date, pbs from gui_project")
    row = c.fetchall();

    if (not row) or (cmd == 'init'):
        c.execute("DELETE FROM gui_project");
        conn.commit();

        # Reading a file
        fid = open("static/data/init_prj.txt", "r");
        for readLine in fid:
            readLine = readLine.strip(' \t\n\r');
            if (len(readLine) > 0):
                tmp = re.split('\t', readLine);
                query = """INSERT INTO gui_project (user_id, name, path1, path2, path3, path4, path5, pbs, date) \
                         VALUES ( "{0}", "{1}", "{2}", "{3}", "{4}", "{5}", "{6}", "{7}", "{8}")""".format( \
                         tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8]);
                #print query;
                c.execute(query);
                conn.commit();
        
        fid.close();
        
        # Retriving updated information
        c.execute("select id, user_id, name, path1, path2, path3, path4, path5, date, pbs from gui_project")
        row = c.fetchall();
        for item in row:
           #print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(item[0], item[1], item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9]);
            parents.append(-1);
            indent.append(0);
            parent.append(-1);
            title.append(item[2]);
            path1.append(item[3]);
            path2.append(item[4]);
            path3.append(item[5]);
            path4.append(item[6]);
            path5.append(item[7]);
            date.append(item[8]);
            rmCrLf = re.sub(r"\s+", "", item[9]);            
            pbs.append(rmCrLf);
            
        for i in range(len(row) ):
            numRec.append(i);
        #print "Number of Items in DB = {}".format(len(row));
        
    elif (cmd == 'delete'):
        #print 'DELETING...';
        c.execute("DELETE FROM gui_project");
        conn.commit();
        
        # display initialized table        
        parents.append(-1);
        indent.append(0);
        parent.append(-1);
        title.append("N/A");
        path1.append("N/A");
        path2.append("N/A");
        path3.append("N/A");
        path4.append("N/A");
        path5.append("N/A");
        pbs.append("N/A");
        date.append("N/A");
        numRec.append(0);           # index starts from 0

    else:
        for item in row:
           #print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(item[1], item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9]);
            parents.append(-1);
            indent.append(0);
            parent.append(-1);
            title.append(item[2]);
            path1.append(item[3]);
            path2.append(item[4]);
            path3.append(item[5]);
            path4.append(item[6]);
            path5.append(item[7]);
            date.append(item[8]);
            
            rmCrLf = re.sub(r"\s+", "", item[9]);
            pbs.append(rmCrLf);

        for i in range(len(row) ):
            numRec.append(i);
        #print "Number of Items in DB = {}".format(len(row));
    

    c = {
            'parents'       : parents,
            'indent'	    : indent,
            'parent'	    : parent,
            'title'	    : title,
            'path1'	    : path1,
            'path2'	    : path2,
            'path3'	    : path3,
            'path4'	    : path4,
            'path5'	    : path5,
            'date'	    : date,
            'pbs'           : pbs,
            'numRec'        : numRec,
        }
    #print c;
    
    conn.close();
    if request.is_ajax() and (request.method == 'POST'):
       return HttpResponse(json.dumps(c));
    else:
        template = 'gui/toyDB_prj.html';
        return render_to_response(template, c, context_instance = RequestContext(request));

def toyView_usr(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
    else:
        cmd = 'http';           # connection with URL

    #---------------------------------------------#
    # HANDLING PROJECT TABLE ** START **
    #---------------------------------------------#
    # Initializing variables
    parents = [];
    indent = [];
    parent = [];
    uid  = [];
    pwd  = [];
    email  = [];
    level  = [];
    numRec = [];
    
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    c.execute("select uid, pwd, email, level from gui_user")
    row = c.fetchall();

    if (not row) or (cmd == 'init'):
        c.execute("DELETE FROM gui_user");
        conn.commit();

        # Reading a file
        fid = open("static/data/init_usr.txt", "r");
        for readLine in fid:
            readLine = readLine.strip(' \t\n\r');
            if (len(readLine) > 0):
                tmp = re.split('\t', readLine);
                #print tmp;
                hpwd = hashlib.md5(tmp[1]);
                hpwd = hpwd.hexdigest();
                query = """INSERT INTO gui_user (uid, pwd, email, level) \
                         VALUES ( "{0}", "{1}", "{2}", {3})""".format( \
                         tmp[0], hpwd, tmp[2], tmp[3]);
                #print query;
                c.execute(query);
                conn.commit();
        
        fid.close();
        
        # Retriving updated information
        c.execute("select uid, pwd, email, level from gui_user")
        row = c.fetchall();
        for item in row:
            parents.append(-1);
            indent.append(0);
            parent.append(-1);
            uid.append(item[0]);
            pwd.append(item[1]);
            email.append(item[2]);
            level.append(item[3]);

        for i in range(len(row) ):
            numRec.append(i);
        #print "Number of Items in DB = {}".format(len(row));
        
    elif (cmd == 'delete'):
        #print 'DELETING...';
        c.execute("DELETE FROM gui_user");
        conn.commit();
        
        # display initialized table        
        parents.append(-1);
        indent.append(0);
        parent.append(-1);
        uid.append("N/A");
        pwd.append("N/A");
        email.append("N/A");
        level.append("N/A");
        numRec.append(0);           # index starts from 0

    else:
        for item in row:
            parents.append(-1);
            indent.append(0);
            parent.append(-1);
            uid.append(item[0]);
            pwd.append(item[1]);
            email.append(item[2]);
            level.append(item[3]);

        for i in range(len(row) ):
            numRec.append(i);
        #print "Number of Items in DB = {}".format(len(row));
    

    c = {
            'parents'       : parents,
            'indent'	    : indent,
            'parent'	    : parent,
            'uid'	    : uid,
            'pwd'	    : pwd,
            'email'	    : email,
            'level'	    : level,
            'numRec'        : numRec,
        }
    #print c;
    
    conn.close();
    if request.is_ajax() and (request.method == 'POST'):
       return HttpResponse(json.dumps(c));
    else:
        template = 'gui/toyDB_usr.html';
        return render_to_response(template, c, context_instance = RequestContext(request));


def prjView(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);        # session will be closed in 10 minutes
    user_id = request.session['user_id'];
    conn = sqlite3.connect(dbName);
    c = conn.cursor();

    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
        pkey  = request.POST.get('pkey');
        pkey = pkey.strip(' \t\n\r');
        if (cmd == 'delete'):
            query = """DELETE FROM gui_project where user_id = "{0}" AND id = {1}""".format(user_id, pkey);
            #print query;
            c.execute(query);
            conn.commit();

    else:
        cmd = 'http';           # connection with URL
        
    #---------------------------------------------#
    # HANDLING PROJECT TABLE ** START **
    #---------------------------------------------#
    # Initializing variables
    pkey = [];
    parents = [];
    indent = [];
    parent = [];
    title  = [];
    path1  = [];
    path2  = [];
    path3  = [];
    path4  = [];
    path5  = [];
    pbs    = [];
    date   = [];
    numRec = [];
    
    c.execute("select id, user_id, name, path1, path2, path3, path4, path5, date, pbs from gui_project")
    row = c.fetchall();

    if (not row):
        # display initialized table
        pkey.append(-1);
        parents.append(-1);
        indent.append(0);
        parent.append(-1);
        title.append("N/A");
        path1.append("N/A");
        path2.append("N/A");
        path3.append("N/A");
        path4.append("N/A");
        path5.append("N/A");
        pbs.append("N/A");
        date.append("N/A");
        numRec.append(0);           # index starts from 0
    else:
        for item in row:
            #print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(item[0], item[1], item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9]);
            #print "DATE:" +  item[8];
            #print "PBS:" + item[9];
            parents.append(-1);
            indent.append(0);
            parent.append(-1);
            pkey.append(item[0]);
            title.append(item[2]);
            path1.append(item[3]);
            path2.append(item[4]);
            path3.append(item[5]);
            path4.append(item[6]);
            path5.append(item[7]);
            date.append(item[8]);
            rmCrLf = re.sub(r"\s+", "", item[9]);
            pbs.append(rmCrLf);
        
        for i in range(len(row) ):
            numRec.append(i);
        #print "Number of Items in DB = {}".format(len(row));
    
    c = {
            'pkey'          : pkey,
            'parents'       : parents,
            'indent'	    : indent,
            'parent'	    : parent,
            'title'	    : title,
            'path1'	    : path1,
            'path2'	    : path2,
            'path3'	    : path3,
            'path4'	    : path4,
            'path5'	    : path5,
            'pbs'           : pbs,
            'date'	    : date,
            'numRec'        : numRec,
        }
    #print c;
    
    conn.close();
    if request.is_ajax() and (request.method == 'POST'):
       return HttpResponse(json.dumps(c));
    else:
        template = 'gui/projectView.html';
        return render_to_response(template, c, context_instance = RequestContext(request));

def prjView_new(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    if request.is_ajax() and (request.method == 'POST'):
        #print "OKAY this is NEW submit!!"
        title    = request.POST.get('title');
        pbs      = request.POST.get('pbs');
        rpaths   = request.POST.getlist('rpaths[]');
        #print rpaths
        out_paths  = request.POST.getlist('out_paths[]');
        #print out_paths
        ex_pythons = request.POST.getlist('ex_pythons[]');
        #print ex_pythons
            
        #pbs.replace("\r", r"\r").replace("\n", r"\n");
        #print pbs
        
        #print request.session.session_key; # get session key
        user = request.session.get('user_id');
        cdate = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
        #cdate = datetime.now().strftime("%A, %d. %B %Y %I:%M%p");
        #print cdate
        
        #fPath = '';
        #vPath = '';
        #print 'learning rpath...'
        #for i in range(len(rpaths)):
        #    fPath = fPath + "path{0}, ".format(i+1);
        #    vPath = vPath + """"{0}", """.format(rpaths[i]);
            
        #query = """INSERT INTO gui_project (user_id, name, {0} date, pbs) \
        #        VALUES ("{1}", "{2}", {3} "{4}", "{5}")""".format(fPath, user, title, vPath, cdate, pbs);
        
        #print dbName
        conn = sqlite3.connect(dbName);
        c = conn.cursor();
        
        # Insert project information into the table
	#print "==== okay insert new value"
        query = """INSERT INTO gui_project (user_id, name, date, pbs) \
                VALUES ("{0}", "{1}", "{2}", "{3}")""".format(user, title, cdate, pbs);
        #print query
        c.execute(query);
        conn.commit();
        
        # Retriving the primary key of project
        query = """SELECT id FROM gui_project WHERE user_id = "{0}" AND date = "{1}" """.format(user, cdate);
        c.execute(query);
        prj_id = c.fetchone();      # primary key = prj_id[0]
        
        #print "PRIMARY KEY IS: {}".format(prj_id[0]);
        
        # Insert input, output, and python path
        for i in range(len(rpaths)):
            query = """INSERT INTO gui_path_input (proj_id, path) \
                       VALUES ({0}, "{1}")""".format(prj_id[0], rpaths[i]);
            #print query;
            c.execute(query);
            conn.commit();

        for i in range(len(out_paths)):
            query = """INSERT INTO gui_path_output (proj_id, path) \
                       VALUES ({0}, "{1}")""".format(prj_id[0], out_paths[i]);
            #print query
            c.execute(query);
            conn.commit();

        for i in range(len(ex_pythons)):
            query = """INSERT INTO gui_path_python (proj_id, path) \
                       VALUES ({0}, "{1}")""".format(prj_id[0], ex_pythons[i]);
            print query
            c.execute(query);
            conn.commit();

        # Display tables;
        
        query = "SELECT id, user_id, name, date, pbs FROM gui_project";
        c.execute(query);
        row = c.fetchall();
        #print query
        #for item in row:
        #    print item;
            
        query = "SELECT id, proj_id, path FROM gui_path_input";
        c.execute(query);
        row = c.fetchall();
        #print query
        #for item in row:
        #    print item;

        query = "SELECT id, proj_id, path FROM gui_path_output";
        c.execute(query);
        row = c.fetchall();
        #print query
        #for item in row:
        #    print item;

        query = "SELECT id, proj_id, path FROM gui_path_python";
        c.execute(query);
        row = c.fetchall();
        #print query
        #for item in row:
        #    print item;
        
        conn.close();
        
    else:
        title       = "";
        pbs         = "";
        rpaths      = "";
        out_paths   = "";
        ex_pythons  = "";
            
        
    c = {
            'title'	    : title,
            'pbs'	    : pbs,
            'rpaths'	    : rpaths,
            'out_paths'     : out_paths,
            'ex_pythons'    : ex_pythons,
        }
        
    if request.is_ajax() and (request.method == 'POST'):
       return HttpResponse(json.dumps(c));
    else:
        template = 'gui/projectView_new.html';
        return render_to_response(template, c, context_instance = RequestContext(request));
    
def prjView_update(request):
    # check out authority
    #print "I'm in prjView_update!"
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);        # session will be closed in 10 minutes
    if request.is_ajax() and (request.method == 'POST'):
        cmd     =  request.POST.get('cmd');
        pkey     = request.POST.get('pkey');
        title    = request.POST.get('title');
        pbs      = request.POST.get('pbs');
        rpaths   = request.POST.getlist('rpaths[]');
        out_paths  = request.POST.getlist('out_paths[]');
        ex_pythons = request.POST.getlist('ex_pythons[]');
        
        user = request.session.get('user_id');
        cdate = datetime.now().strftime("%Y-%m-%d %H:%M:%S");

        if cmd == 'reload':            
            #print dbName
            conn = sqlite3.connect(dbName);
            c = conn.cursor();
            
            # Update project information
            query = """UPDATE gui_project SET name="{0}", date="{1}", pbs="{2}" \
                    WHERE id={3} AND user_id="{4}" """.format(title, cdate, pbs, pkey, user);
            #print query
            c.execute(query);
            conn.commit();
            
            # Deleting associated tables - this makes easier to update it
            query = "DELETE FROM gui_path_input WHERE proj_id={0}".format(pkey);
            #print query
            c.execute(query);
            conn.commit();
            
            # Insert input, output, and python path
            for i in range(len(rpaths)):
                query = """INSERT INTO gui_path_input (proj_id, path) \
                           VALUES ({0}, "{1}")""".format(pkey, rpaths[i]);
                #print query
                c.execute(query);
                conn.commit();
    
            # Deleting associated tables - this makes easier to update it
            query = "DELETE FROM gui_path_output WHERE proj_id={0}".format(pkey);
            #print query
            c.execute(query);
            conn.commit();
    
            for i in range(len(out_paths)):
                query = """INSERT INTO gui_path_output (proj_id, path) \
                           VALUES ({0}, "{1}")""".format(pkey, out_paths[i]);
                #print query
                c.execute(query);
                conn.commit();
    
            # Deleting associated tables - this makes easier to update it
            query = "DELETE FROM gui_path_python WHERE proj_id={0}".format(pkey);
            #print query
            c.execute(query);
            conn.commit();
    
            for i in range(len(ex_pythons)):
                query = """INSERT INTO gui_path_python (proj_id, path) \
                           VALUES ({0}, "{1}")""".format(pkey, ex_pythons[i]);
                #print query
                c.execute(query);
                conn.commit();
    else:
        pkey        = "";
        title       = "";
        pbs         = "";
        rpaths      = "";
        out_paths   = "";
        ex_pythons  = "";
            
        
    c = {
            'pkey'          : pkey,
            'title'	    : title,
            'pbs'	    : pbs,
            'rpaths'	    : rpaths,
            'out_paths'     : out_paths,
            'ex_pythons'    : ex_pythons,
        }
        
    if request.is_ajax() and (request.method == 'POST'):
       return HttpResponse(json.dumps(c));
    else:
        template = 'gui/projectView_new.html';
        return render_to_response(template, c, context_instance = RequestContext(request));
        

def prjView_delete(request):
    # check out authority
    #print "I'm in prjView_delete!"
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);        # session will be closed in 10 minutes
    if request.is_ajax() and (request.method == 'POST'):
        cmd     =  request.POST.get('cmd');
        pkey     = request.POST.get('pkey');
        
        user = request.session.get('user_id');
        cdate = datetime.now().strftime("%Y-%m-%d %H:%M:%S");

        if cmd == 'prjDelete':            
            #print dbName
            conn = sqlite3.connect(dbName);
            c = conn.cursor();
            
            # DELETING project information
            query = "DELETE FROM gui_project WHERE id={0}".format(pkey);
            #print query
            c.execute(query);
            conn.commit();

            query = "DELETE FROM gui_job WHERE proj_id={0}".format(pkey);
            #print query
            c.execute(query);
            conn.commit();

            query = "DELETE FROM gui_path_input WHERE proj_id={0}".format(pkey);
            #print query
            c.execute(query);
            conn.commit();

            query = "DELETE FROM gui_path_output WHERE proj_id={0}".format(pkey);
            #print query
            c.execute(query);
            conn.commit();

            query = "DELETE FROM gui_path_python WHERE proj_id={0}".format(pkey);
            #print query
            c.execute(query);
            conn.commit();
            
            conn.close();
    else:
        pkey        = "";
            
        
    c = {
            'pkey'          : pkey,
        }
        
    if request.is_ajax() and (request.method == 'POST'):
       return HttpResponse(json.dumps(c));
    else:
        template = 'gui/projectView_new.html';
        return render_to_response(template, c, context_instance = RequestContext(request));
        
    

def jobView(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )
    
    request.session.set_expiry(SESSION_TIME_OUT);
    template = 'gui/jobView.html';
    c = {
        'cmd': 'init',
    }
    return render_to_response(template, c, context_instance = RequestContext(request));

def jobView_jqGrid_para(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )
    
    request.session.set_expiry(SESSION_TIME_OUT);
    job_id = request.GET.get('job_id');

    if request.GET.get("sidx"):
        sidx = request.GET["sidx"];
        sord = request.GET["sord"];
        cpage   = request.GET["page"];
        cpage   = int(cpage);
        maxrow  = request.GET["rows"];
        maxrow  = int(maxrow);    
        
    if cpage == 1:
        st_index = 0;
        ed_index = maxrow;
    else:
        st_index = maxrow * cpage - maxrow;
        ed_index = maxrow * cpage;
        
        
    if (request.GET["_search"] == 'true'):
        sfilter = request.GET["filters"];
        search_dic =json.loads(sfilter);      # convert string to dictionary
        search_rules = search_dic["rules"][0];
        search_field = search_rules["field"];
        search_data = search_rules["data"];
        query = """select id, anaz, para, val, status from gui_parameter where job_id = {0} AND {1} like "%{2}%" order by {3} {4} limit {5}, {6}""".format(job_id, search_field, search_data, sidx, sord, st_index, ed_index);
        #print query
    else:   
        query = "select id, anaz, para, val, status from gui_parameter where job_id = {0} order by {1} {2}".format(job_id, sidx, sord);
        #print query

    #query = "select id, job_id, anaz, para, val, status from gui_parameter where job_id = {0} order by {1} {2}".format(job_id, sidx, sord);
    #query = "select id, anaz, para, val, status from gui_parameter where job_id = {0} order by {1} {2}".format(job_id, sidx, sord);
    #print query;
    conn = sqlite3.connect(dbName);
    c = conn.cursor();    
    c.execute(query);
    row = c.fetchall();
    conn.close();
    
    rows = [];
    row_cnt = 0;
    if st_index == 0:    
	for item in row:
	    row_cnt = row_cnt + 1;
	    tmp = { "id": str(row_cnt), "cell": item}
	    rows.append(tmp);
    else:
	for item in row:
	    row_cnt = row_cnt + 1;
	    if (row_cnt > st_index):
		tmp = { "id": str(row_cnt), "cell": item}
		rows.append(tmp);

    # calculating total number of pages
    tpages = int(math.ceil(float(row_cnt) / float(maxrow)));
    #print "tpages: {}".format(tpages);
    if tpages < 1:
        tpages = 1;
    
    total_pages = tpages;

    c = {
            'total'	: total_pages,
            'page'      : cpage,
            'records'   : row_cnt,
            'rows'      : rows,            
        }
    return HttpResponse(json.dumps(c));

    
def jobView_jqGrid_job(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )
    
    request.session.set_expiry(SESSION_TIME_OUT);
    prj_id = request.GET.get('prj_id');
    
    if request.GET.get("sidx"):
        sidx = request.GET["sidx"];
        sord = request.GET["sord"];
        cpage   = request.GET["page"];
        cpage   = int(cpage);
        maxrow  = request.GET["rows"];
        maxrow  = int(maxrow);    
        
    if cpage == 1:
        st_index = 0;
        ed_index = maxrow;
    else:
        st_index = maxrow * cpage - maxrow;
        ed_index = maxrow * cpage;

    if (request.GET["_search"] == 'true'):
        sfilter = request.GET["filters"];
        search_dic =json.loads(sfilter);      # convert string to dictionary
        search_rules = search_dic["rules"][0];
        search_field = search_rules["field"];
        search_data = search_rules["data"];
        query = """select id, name, anaz, status, output, stime, etime from gui_job where proj_id = {0} AND {1} like "%{2}%" order by {3} {4} limit {5}, {6}""".format(prj_id, search_field, search_data, sidx, sord, st_index, ed_index);
        #print query
    else:   
        query = "select id, name, anaz, status, output, stime, etime from gui_job where proj_id = {0} order by {1} {2}".format(prj_id, sidx, sord);
        #print query

    #query = "select id, proj_id, name, anaz, status, output, stime, etime from gui_job where proj_id = {0} order by {1} {2}".format(prj_id, sidx, sord);
    #print query;
    conn = sqlite3.connect(dbName);
    c = conn.cursor();    
    c.execute(query);
    row = c.fetchall();
    conn.close();
    
    rows = [];
    row_cnt = 0;
    
    if st_index == 0:
	for item in row:
	    row_cnt = row_cnt + 1;
	    tmp = { "id": str(row_cnt), "cell": item}
	    rows.append(tmp);
    else:
	for item in row:
	    row_cnt = row_cnt + 1;
	    if (row_cnt > st_index):
		tmp = { "id": str(row_cnt), "cell": item}
		rows.append(tmp);
	
        
    # calculating total number of pages
    tpages = int(math.ceil(float(row_cnt) / float(maxrow)));
    #print "tpages: {}".format(tpages);
    if tpages < 1:
        tpages = 1;
    
    total_pages = tpages;

    c = {
            'total'	: total_pages,
            'page'      : cpage,
            'records'   : row_cnt,
            'rows'      : rows,            
        }
    return HttpResponse(json.dumps(c));

    

def jobView_jqGrid_prj(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    user_id = request.session['user_id'];
    #print request.GET
        
    if request.GET.get("sidx"):
        sidx    = request.GET["sidx"];
        sord    = request.GET["sord"];
        cpage   = request.GET["page"];
        cpage   = int(cpage);
        maxrow  = request.GET["rows"];
        maxrow  = int(maxrow);    
        
    if cpage == 1:
        st_index = 0;
        ed_index = maxrow;
    else:
        st_index = maxrow * cpage - maxrow;
        ed_index = maxrow * cpage;
    
    #if request.is_ajax() and (request.method == 'POST'):
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    
    if (request.GET["_search"] == 'true'):
        sfilter = request.GET["filters"];
        search_dic =json.loads(sfilter);      # convert string to dictionary
        search_rules = search_dic["rules"][0];
        search_field = search_rules["field"];
        search_data = search_rules["data"];
        query = """select id, user_id, name, date, pbs from gui_project where user_id = "{0}" AND {1} like "%{2}%" order by {3} {4} limit {5}, {6}""".format(user_id, search_field, search_data, sidx, sord, st_index, ed_index);
        #print query
    else:   
        query = "select id, user_id, name, date, pbs from gui_project where user_id = '{0}' order by {1} {2}".format(user_id, sidx, sord);
        #print query
    
    c.execute(query);
    row = c.fetchall();
    conn.close();

    rows = [];
    row_cnt = 0;
    if st_index == 0:
	for item in row:
	    row_cnt = row_cnt + 1;
	    tmp = { "id": str(row_cnt), "cell": item}
	    rows.append(tmp);
    else:
	for item in row:
	    row_cnt = row_cnt + 1;
	    if (row_cnt > st_index):
		tmp = { "id": str(row_cnt), "cell": item}
		rows.append(tmp);
    
    # calculating total number of pages
    tpages = int(math.ceil(float(row_cnt) / float(maxrow)));
    #print "tpages: {}".format(tpages);
    if tpages < 1:
        tpages = 1;
    
    total_pages = tpages;

    c = {
            'total'	: total_pages,
            'page'      : cpage,
            'records'   : row_cnt,
            'rows'      : rows,            
        }
    return HttpResponse(json.dumps(c));

def resultView_jqGrid_results(request):
    #print "OKAY i am In RESUT VIEW"
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )
    
    request.session.set_expiry(SESSION_TIME_OUT);
    prj_id = request.GET.get('prj_id');
    
    if request.GET.get("sidx"):
        sidx = request.GET["sidx"];
        sord = request.GET["sord"];
        cpage   = request.GET["page"];
        cpage   = int(cpage);
        maxrow  = request.GET["rows"];
        maxrow  = int(maxrow);    
        
    if cpage == 1:
        st_index = 0;
        ed_index = maxrow;
    else:
        st_index = maxrow * cpage - maxrow;
        ed_index = maxrow * cpage;

    if (request.GET["_search"] == 'true'):
	#print "--- search true at resultView ----"
        sfilter = request.GET["filters"];
        search_dic =json.loads(sfilter);      # convert string to dictionary
        search_rules = search_dic["rules"][0];
        search_field = search_rules["field"];
        search_data = search_rules["data"];
	# select all job ID corresponding to 'prj_id'
        query = """select id from gui_job where proj_id = {0} order by {1} {2}""".format(prj_id, sidx, sord);
	#print query
	conn = sqlite3.connect(dbName);
	c = conn.cursor();    
	c.execute(query);
	tmp = c.fetchall();
	JOBs = [];
	for job_id in tmp:
	    query = """select id, job_id, name, img, txt, gzip from gui_outputs where job_id = {0} AND {1} like "%{2}%" order by {3} {4}""".format(job_id[0], search_field, search_data, sidx, sord);
	    #print query
	    c.execute(query);
	    row = c.fetchall();
	    JOBs.append(row);
	conn.close();
	#print JOBs
	conn.close();
	
    else:
	#print "--- search FALSE at resultView ----"
	JOBs = [];
	query = """select id from gui_job where proj_id = {0} order by {1} {2}""".format(prj_id, sidx, sord);
	#print query
	conn = sqlite3.connect(dbName);
	c = conn.cursor();    
	c.execute(query);
	tmp = c.fetchall();	# retrieve all related jobs corresponding to the current project 
	for job_id in tmp:
	    query = """select id, job_id, name, img, txt, gzip from gui_outputs where job_id = {0} order by {1} {2}""".format(job_id[0], sidx, sord);
	    #print query
	    c.execute(query);
	    row = c.fetchall();
	    #print row
	    JOBs.append(row);
	conn.close();
	#print row
	#print "============== UPPER ROW, BELOW JOBs ===============";
	#print JOBs

    rows = [];
    row_cnt = 0;
    if st_index == 0:
	for row in JOBs:
	    for item in row:
		#item_arr = [ i for i in item ];
		row_cnt = row_cnt + 1;
		tmp = { "id": str(row_cnt), "cell": item}
		#print tmp;
		rows.append(tmp);
    else:
	for row in JOBs:
	    for item in row:
		row_cnt = row_cnt + 1;
		if (row_cnt > st_index):
		    tmp = { "id": str(row_cnt), "cell": item}
		    rows.append(tmp);
    
        
    # calculating total number of pages
    tpages = int(math.ceil(float(row_cnt) / float(maxrow)));
    #print "tpages: {}".format(tpages);
    if tpages < 1:
        tpages = 1;
    
    total_pages = tpages;

    c = {
            'total'	: total_pages,
            'page'      : cpage,
            'records'   : row_cnt,
            'rows'      : rows,            
        }
    return HttpResponse(json.dumps(c));

    

def resultView_jqGrid_del_results(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    user_id = request.session['user_id'];
    
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
    else:
        cmd = 'http';           # connection with URL
        c = {
                    'errMsg'	    : 'Session has been expired!',
            }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request))
    
    if cmd == 'del_user':
        #print "OKAY I am in del_user"
        # Extract my level
        conn = sqlite3.connect(dbName);
        c = conn.cursor();
        query = "select level from gui_user where uid='{0}'".format(user_id);
        c.execute(query);
        my_level = c.fetchone();
        my_level = int(my_level[0]);
        delUsers = request.POST.getlist('uid[]');
        fUsers = [];                                 # the list of users failed for deletion 
        for t_user in delUsers:
            query = "select level from gui_user where uid='{0}'".format(t_user);
            c.execute(query);
            user_level = c.fetchone();
            user_level = int(user_level[0]);
            if (my_level > user_level) or (user_id == t_user):
                query = "delete from gui_user where uid='{}'".format(t_user);
                #print query
                c.execute(query);
                conn.commit();
            else:
                #print t_user
                fUsers.append(t_user);
        conn.close();
        c = {
            'fUsers'    : fUsers,
        }
    return HttpResponse(json.dumps(c));

    
def userView(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
    else:
        cmd = 'http';           # connection with URL

    #---------------------------------------------#
    # HANDLING PROJECT TABLE ** START **
    #---------------------------------------------#
    # Initializing variables
    parents = [];
    indent = [];
    parent = [];
    uid  = [];
    pwd  = [];
    email  = [];
    level  = [];
    numRec = [];
    
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    c.execute("select uid, pwd, email, level from gui_user")
    row = c.fetchall();

    if (not row):
        # display initialized table        
        parents.append(-1);
        indent.append(0);
        parent.append(-1);
        uid.append("N/A");
        pwd.append("N/A");
        email.append("N/A");
        level.append("N/A");
        numRec.append(0);           # index starts from 0
    else:
        for item in row:
            parents.append(-1);
            indent.append(0);
            parent.append(-1);
            uid.append(item[0]);
            pwd.append(item[1]);
            email.append(item[2]);
            level.append(item[3]);

        for i in range(len(row) ):
            numRec.append(i);
        #print "Number of Items in DB = {}".format(len(row));
    

    c = {
            'parents'       : parents,
            'indent'	    : indent,
            'parent'	    : parent,
            'uid'	    : uid,
            'pwd'	    : pwd,
            'email'	    : email,
            'level'	    : level,
            'numRec'        : numRec,
        }
    #print c;
    
    conn.close();
    if request.is_ajax() and (request.method == 'POST'):
       return HttpResponse(json.dumps(c));
    else:
        template = 'gui/user.html';
        return render_to_response(template, c, context_instance = RequestContext(request));

def usrView_jqGrid_create_user(request):
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    user_id = request.session['user_id'];
    
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
    else:
        cmd = 'http';           # connection with URL
        c = {
                    'errMsg'	    : 'Session has been expired!',
            }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request))
    if cmd == 'get_level':
        level = getUsrLevel(dbName, user_id);
        c = {
            'user' : user_id,
            'level': level,
        }
    
    if cmd == 'newUser':
        new_user    = request.POST.get('user');
        new_pwd     = request.POST.get('pwd');
        hpwd = hashlib.md5(new_pwd);
        hpwd = hpwd.hexdigest();
        new_email   = request.POST.get('email');
        new_level   = request.POST.get('level');
        
        conn = sqlite3.connect(dbName);
        c = conn.cursor();
        query = "INSERT INTO gui_user (uid, pwd, email, level) VALUES ('{0}', '{1}', '{2}', {3})".format(new_user, hpwd, new_email, new_level)
        c.execute(query);
        conn.commit();
        conn.close();
        c = {
            'user'	: new_user,
            'pwd'	: hpwd,
            'email'	: new_email,
            'level'	: new_level,
        }
        
    if cmd == 'editUser':
        #print "I am in EditUSER"
        edit_user    = request.POST.get('user');
        edit_pwd     = request.POST.get('pwd');
        hpwd = hashlib.md5(edit_pwd);
        hpwd = hpwd.hexdigest();
        edit_email   = request.POST.get('email');
        edit_level   = request.POST.get('level');

        conn = sqlite3.connect(dbName);
        c = conn.cursor();
        query = "UPDATE gui_user SET uid='{0}', pwd='{1}', email='{2}', level={3} WHERE uid='{4}' ".format(edit_user, hpwd, edit_email, edit_level, edit_user);
        #print query
        c.execute(query);
        conn.commit();
        conn.close();
        c = {
            'user'	: edit_user,
            'pwd'	: hpwd,
            'email'	: edit_email,
            'level'	: edit_level,
        }
        
    return HttpResponse(json.dumps(c));


def usrView_jqGrid_del_user(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    user_id = request.session['user_id'];
    
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
    else:
        cmd = 'http';           # connection with URL
        c = {
                    'errMsg'	    : 'Session has been expired!',
            }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request))
    
    if cmd == 'del_user':
        #print "OKAY I am in del_user"
        # Extract my level
        conn = sqlite3.connect(dbName);
        c = conn.cursor();
        query = "select level from gui_user where uid='{0}'".format(user_id);
        c.execute(query);
        my_level = c.fetchone();
        my_level = int(my_level[0]);
        delUsers = request.POST.getlist('uid[]');
        fUsers = [];                                 # the list of users failed for deletion 
        for t_user in delUsers:
            query = "select level from gui_user where uid='{0}'".format(t_user);
            c.execute(query);
            user_level = c.fetchone();
            user_level = int(user_level[0]);
            if (my_level > user_level) or (user_id == t_user):
                query = "delete from gui_user where uid='{}'".format(t_user);
                #print query
                c.execute(query);
                conn.commit();
            else:
                #print t_user
                fUsers.append(t_user);
        conn.close();
        c = {
            'fUsers'    : fUsers,
        }
    return HttpResponse(json.dumps(c));



def usrView_jqGrid_usr(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    user_id = request.session['user_id'];
    #print request.GET
        
    if request.GET.get("sidx"):
        sidx    = request.GET["sidx"];
        sord    = request.GET["sord"];
        cpage   = request.GET["page"];
        cpage   = int(cpage);
        maxrow  = request.GET["rows"];
        maxrow  = int(maxrow);    
        
    if cpage == 1:
        st_index = 0;
        ed_index = maxrow;
    else:
        st_index = maxrow * cpage - maxrow;
        ed_index = maxrow * cpage;
    
    #if request.is_ajax() and (request.method == 'POST'):
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    
    if (request.GET["_search"] == 'true'):
        sfilter = request.GET["filters"];
        search_dic =json.loads(sfilter);      # convert string to dictionary
        search_rules = search_dic["rules"][0];
        search_field = search_rules["field"];
        search_data = search_rules["data"];
        query = """select uid, email, level from gui_user where {0} like "%{1}%" order by {2} {3}""".format(search_field, search_data, sidx, sord);
        #print query
    else:   
        query = "select uid, email, level from gui_user order by {0} {1}".format(sidx, sord);
        #print query
    
    c.execute(query);
    row = c.fetchall();
    conn.close();

    rows = [];
    row_cnt = 0;
    if st_index == 0:
	for item in row:
	    row_cnt = row_cnt + 1;
	    tmp = { "id": str(row_cnt), "cell": item}
	    rows.append(tmp);
    else:
	for item in row:
	    row_cnt = row_cnt + 1;
	    if (row_cnt > st_index):
		tmp = { "id": str(row_cnt), "cell": item}
		rows.append(tmp);
	
    
    # calculating total number of pages
    tpages = int(math.ceil(float(row_cnt) / float(maxrow)));
    #print "tpages: {}".format(tpages);
    if tpages < 1:
        tpages = 1;
    
    total_pages = tpages;
    #print rows
    c = {
            'total'	: total_pages,
            'page'      : cpage,
            'records'   : row_cnt,
            'rows'      : rows,            
        }
    return HttpResponse(json.dumps(c));

def stanalyzer(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);        # session will be closed in 10 minutes
    user_id = request.session['user_id'];
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    #print "OKAY~~~~~~~~ I am here !!!!!!!!!!!!!!!!!!!!"
    if request.is_ajax() and (request.method == 'POST'):
        #print "This is AJAX!"
        cmd   = request.POST.get('cmd');
        #print "Okay you requested CMD as {}".format(cmd);
        if (cmd == 'path'):
            pkey  = request.POST.get('pkey');
            pkey = pkey.strip(' \t\n\r');

            #query = "select id, user_id, name, path1, path2, path3, path4, path5, date, pbs from gui_project where id = {0}".format(pkey);
            query = "select id, user_id, name, date, pbs from gui_project where user_id = '{0}' AND id = {1}".format(user_id, pkey);
            c.execute(query);
            prj = c.fetchall();
            
            query = "select id, proj_id, path from gui_path_input where proj_id = {0}".format(pkey);
            c.execute(query);
            path_in = c.fetchall();
            
            query = "select id, proj_id, path from gui_path_output where proj_id = {0}".format(pkey);
            c.execute(query);
            path_out = c.fetchall();

            query = "select id, proj_id, path from gui_path_python where proj_id = {0}".format(pkey);
            c.execute(query);
            path_py = c.fetchall();
            
        elif (cmd =='reload'):
            #print "This is  AJAX with relaod!!!!!!"
            #c.execute("select id, user_id, name, path1, path2, path3, path4, path5, date, pbs from gui_project order by date desc")
            # find most recent one
            query = "select id, user_id, name, date, pbs from gui_project where user_id = '{0}' order by date desc".format(user_id);
            c.execute(query);
            prj = c.fetchall();
            pkey = prj[0][0];
            
            query = "select id, proj_id, path from gui_path_input where proj_id = {0}".format(pkey);
            c.execute(query);
            path_in = c.fetchall();
            
            query = "select id, proj_id, path from gui_path_output where proj_id = {0}".format(pkey);
            c.execute(query);
            path_out = c.fetchall();
    
            query = "select id, proj_id, path from gui_path_python where proj_id = {0}".format(pkey);
            c.execute(query);
            path_py = c.fetchall();
        elif (cmd == 'prjDelete'):
            #print "CMD: prjDelte"
            query = "select id, user_id, name, date, pbs from gui_project where user_id = '{0}' order by date desc".format(user_id);
            c.execute(query);
            prj = c.fetchall();
            if len(prj) < 1 :
                pkey = ['-1'];
                path_in  = ['N/A'];
                path_out = ['N/A'];
                path_py  = ['N/A'];
            else:
                pkey = prj[0][0];
            
                query = "select id, proj_id, path from gui_path_input where proj_id = {0}".format(pkey);
                c.execute(query);
                path_in = c.fetchall();
                
                query = "select id, proj_id, path from gui_path_output where proj_id = {0}".format(pkey);
                c.execute(query);
                path_out = c.fetchall();
        
                query = "select id, proj_id, path from gui_path_python where proj_id = {0}".format(pkey);
                c.execute(query);
                path_py = c.fetchall();
            
    else:
        #print "This is NOT AJAX!!!!!!"
        #c.execute("select id, user_id, name, path1, path2, path3, path4, path5, date, pbs from gui_project order by date desc")
        # find most recent one
        query = "select id, user_id, name, date, pbs from gui_project where user_id = '{0}' order by date desc".format(user_id);
        c.execute(query);
        prj = c.fetchall();
        if len(prj) < 1 :
            pkey = ['-1'];
            path_in  = ['N/A'];
            path_out = ['N/A'];
            path_py  = ['N/A'];
        else:
            pkey = prj[0][0];
        
            query = "select id, proj_id, path from gui_path_input where proj_id = {0}".format(pkey);
            c.execute(query);
            path_in = c.fetchall();
            
            query = "select id, proj_id, path from gui_path_output where proj_id = {0}".format(pkey);
            c.execute(query);
            path_out = c.fetchall();
    
            query = "select id, proj_id, path from gui_path_python where proj_id = {0}".format(pkey);
            c.execute(query);
            path_py = c.fetchall();

    #print "length ROW: {}".format(len(row));
    #for item in row:
    #    print "id:{0}\nuser_id:{1}\nTitle:{2}\npath1:{3}\npath2:{4}\npath3:{5}\npath4:{6}\npath5:{7}\ndate:{8}\npbs:{9}\n".format(item[0], item[1], item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9]);

    # Initializing variables
    pkey = [];
    user = [];
    title  = [];
    ptitle= zip(pkey, title);
    
    path_inputs_id = [];
    path_inputs_proj_id = [];
    path_inputs_path = [];
    
    path_outputs_id = [];
    path_outputs_proj_id = [];
    path_outputs_path = [];

    path_pythons_id = [];
    path_pythons_proj_id = [];
    path_pythons_path = [];
    
    pbs    = [];
    date   = [];
    numRec = [];
    fList  = [];    # file list in the base directory
    if (not prj):
        # display initialized table
        pkey.append("");
        user.append("");
        title.append("");
        path_inputs_id.append("");
        path_inputs_proj_id.append("");
        path_inputs_path.append("");
        path_outputs_id.append("");
        path_outputs_proj_id.append("");
        path_outputs_path.append("");
        path_pythons_id.append("");
        path_pythons_proj_id.append("");
        path_pythons_path.append("");
        pbs.append("");
        date.append("");
        numRec.append("");           # index starts from 0
        fList.append("");
    else:
        # print "ROW = {0}".format(len(row));
        # for project
        for item in prj:
            pkey.append(item[0]);
            user.append(item[1]);
            title.append(item[2]);
            date.append(item[3]);
            pbs.append(item[4]);
            
        for item in path_in:
            path_inputs_id.append(item[0]);
            path_inputs_proj_id.append(item[1]);
            path_inputs_path.append(item[2]);
            
        for item in path_out:
            path_outputs_id.append(item[0]);
            path_outputs_proj_id.append(item[1]);
            path_outputs_path.append(item[2]);
            
        for item in path_py:
            path_pythons_id.append(item[0]);
            path_pythons_proj_id.append(item[1]);
            path_pythons_path.append(item[2]);

        #print "Number of Items in DB = {}".format(len(row));
        for i in range(len(prj)):
            numRec.append(i);

        server = serverside(path_inputs_path[0]);
        fList = server.showdir();
        """
        status =  server.isvalidate();
        #print "Status for {} = {}".format(path1[0], status);
        if not status:
            fList = server.showdir();
        else:
            fList = [];
        """ 
        ptitle = zip(pkey, title);
    c = {
            'PROJECT_ROOT'          : PROJECT_ROOT,
            'dbName'                : dbName,
            'MEDIA_HOME'            : MEDIA_HOME,
            'ANALYZER_HOME'         : ANALYZER_HOME,
            'user'                  : request.session['user_id'],
            'pkey'                  : pkey,
            'title'	            : title,
            'ptitle'                : ptitle,
            'path_inputs_id'        : path_inputs_id,
            'path_inputs_proj_id'   : path_inputs_proj_id,
            'path_inputs_path'      : path_inputs_path,
            'path_outputs_id'       : path_outputs_id,
            'path_outputs_proj_id'  : path_outputs_proj_id,
            'path_outputs_path'     : path_outputs_path,
            'path_pythons_id'       : path_pythons_id,
            'path_pythons_proj_id'  : path_pythons_proj_id,
            'path_pythons_path'     : path_pythons_path,
            'pbs'                   : pbs,
            'date'	            : date,
            'numRec'                : numRec,
            'fList'                 : fList,
        }

 
    conn.close();

    if request.is_ajax() and (request.method == 'POST'):
        return HttpResponse(json.dumps(c));
    else:
        #template = 'gui/stanalyzer.html';
        template = 'desktop/index.html';
        return render_to_response(template, c, context_instance = RequestContext(request));



def stanalyzer_info(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);        # session will be closed in 10 minutes
    fList = 'File not found!';
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
        
        # ----------------------------------------
        # Get list of files based on the directory
        # ----------------------------------------
        if (cmd == 'get_flist'):
            # get file list for the given path
            path = request.POST.get('path');
            path = eval_path(path);
            #if (path[len(path)-1] != '/'):
            #    path = path + '/';
            #else:
                #print "passed!";
                #print "path={}, End with [{}]".format(path, path[len(path)-1]);
            server = serverside(path);
            fList = server.showdir();
            c = {
                    'fList' : fList,
                }

        # ---------------------------------------------
        # Processing information requried for job submitting
        # ---------------------------------------------
        if (cmd == 'get_structure'):
            pkey    = request.POST.get('pkey');
            title   = request.POST.get('title');
            bpath   = request.POST.get('bpath');
            stfile  = request.POST.get('stfile');
            trjFile = request.POST.getlist('trjFile[]');
    
            # --- removing white spaces ---
            pkey = pkey.strip(' \t\n\r');
            title = title.strip(' \t\n\r');
            bpath = bpath.strip(' \t\n\r');
            stfile = stfile.strip(' \t\n\r');

            psf = '{0}/{1}'.format(bpath, stfile);
            trj = '{0}/{1}'.format(bpath, trjFile[0]);
            
            #print psf
            #print trj
            # call MDsimulation 
            objSim = simulation(psf, trj);
            
            # get frame, atom, and segment information            
            num_frm = objSim.num_frm;
            num_atom = objSim.num_atom;
            segList = objSim.segList;
            
            # assume the first segments are chosen
            CresList = objSim.get_seg_residues(segList[0]);
            #print '--- CresList ---'
            #print CresList
            
            resList = [];
            resID   = [];
            for i in range(len(CresList)):
                resList.append(CresList.residues[i].name);
                resID.append(CresList.residues[i].id);
            
            #print '------resList-------'
            #print resList
            
            c = {
                'bpath'	    : bpath,
                'trjFile'   : trjFile,
                'stfile'    : stfile,
                'num_frm'   : num_frm,
                'num_atom'  : num_atom,
                'segList'   : segList,
                'resList'   : resList,
                'resID'	    : resID,
            }
            

        # ---------------------------------------------
        # Get segment information
        # ---------------------------------------------
        if (cmd == 'get_segment'):
            segid      = request.POST.get('segID');
            bpath      = request.POST.get('bpath');
            stfile     = request.POST.get('stfile');
            trjFile    = request.POST.getlist('trjFile[]');

            # --- removing white spaces ---
            segid  = segid.strip(' \t\n\r');
            bpath  = bpath.strip(' \t\n\r');
            stfile = stfile.strip(' \t\n\r');
            
            # ---- checkout path: make sure end with '/'
            bpath = eval_path(bpath);

            #print segid
            #print bpath
            #print stfile
            #print trjFile[0]
            
            psf = '{0}/{1}'.format(bpath, stfile);
            trj = '{0}/{1}'.format(bpath, trjFile[0]);
            #print psf
            #print trj
            # call MDanalysis
            objSim = simulation(psf, trj);
            #print 'Object created'
            seg_name = objSim.get_segname(segid);
            #print seg_name
            CresInfo   = objSim.get_seg_residues(seg_name);
            resList = [];
            resID =[];
            for i in range(len(CresInfo)):
                resList.append(CresInfo[i].name);
                resID.append(CresInfo[i].id);
                
            #print '---resList---'
            #print resList
            #print '---resID---'
            #print resID
            
            c = {
                'resList' : resList,
                'resID'   : resID,
            }

        return HttpResponse(json.dumps(c));
    else:
        return HttpResponseRedirect("/gui/")





def stanalyzer_sendJob(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    data = request.POST;
    
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
        
        # ----------------------------------------
        # Writing Job script and send it
        # ----------------------------------------
        if (cmd == 'sendJob'):
            job_title   = request.POST.get('job_title');
            pkey        = request.POST.get('pkey');
            ptitle      = request.POST.get('ptitle');
            machine     = request.POST.get('machine');
            func_name   = request.POST.getlist('funcName[]');
            Paras       = request.POST.getlist('Paras[]');
            ParaInfo    = request.POST.getlist('ParaInfo[]');
            bpath       = request.POST.get('bpath');
            stfile      = request.POST.get('stfile');
            path_output = request.POST.get('path_output');
            path_python = request.POST.get('path_python');
            trjFile     = request.POST.getlist('trjFile[]');
            pbs         = request.POST.get('pbs');

            # --- removing white spaces ---
            bpath  = bpath.strip(' \t\n\r');
            stfile = stfile.strip(' \t\n\r');
            path_output = path_output.strip(' \t\n\r');
            path_python = path_python.strip(' \t\n\r');
            
            # ---- checkout path: make sure end with '/'
            bpath = eval_path(bpath);
            path_output = eval_path(path_output);
            
            #--- parsing List ----
            Paras = parseWrapList(Paras);
            ParaInfo = parseWrapList(ParaInfo);
            
            #PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__));

            #--- Get unique session ID
            user = request.session['user_id'];
            date = datetime.now().strftime("%Y%m%d%H%M%S%f");
            rndStr = rand_N_letters(6);
            
            #-------------------------------------------
            # creating SESSION_HOME directory
            #-------------------------------------------
            #print "creating directory"
            #print __file__
            
            SESSION_HOME = "{}/media/{}".format(PROJECT_ROOT[0:len(PROJECT_ROOT)-4],user);
            OUTPUT_HOME  = eval_path(SESSION_HOME);
            if (path_output == OUTPUT_HOME):
                #print SESSION_HOME
                if not (os.path.isdir(SESSION_HOME)):
                    #print "Creating directory into {}".format(SESSION_HOME)
                    os.mkdir(SESSION_HOME);
                    
                SESSION_HOME = "{}/{}{}".format(SESSION_HOME, date, rndStr);
                if not (os.path.isdir(SESSION_HOME)):
                    #print "Creating directory into {}".format(SESSION_HOME)
                    os.mkdir(SESSION_HOME);
                else:
                    while not (os.path.isdir(SESSION_HOME)):
                        rndStr = rand_N_letters(6);
                        SESSION_HOME = "{}/{}{}".format(SESSION_HOME, date, rndStr);
                        if not (os.path.isdir(SESSION_HOME)):
                            #print "Creating directory into {}".format(SESSION_HOME)
                            os.mkdir(SESSION_HOME);
                            break;
                OUTPUT_HOME = eval_path(SESSION_HOME);
            else:
                OUTPUT_HOME = path_output;
                
            PBS_HOME = "{}pbs".format(OUTPUT_HOME);
            if not (os.path.isdir(PBS_HOME)):
                #print "creating PBS_HOME: {}".format(PBS_HOME);
                os.mkdir(PBS_HOME);

            SH_HOME  = "{}sh".format(OUTPUT_HOME);
            if not (os.path.isdir(SH_HOME)):
                #print "creating SH_HOME: {}".format(SH_HOME);
                os.mkdir(SH_HOME);
            
            # Wriring variables into binary format
            file_out = "{}para".format(OUTPUT_HOME);
            
            # ---------------------------- Writing PBS
            #print "PBS: {}".format(pbs);
            for ifunc in func_name:
                tmp = "{}/{}.pbs".format(PBS_HOME, ifunc);
                fid_i = open(tmp, 'w');
                fid_i.write(pbs);
                
                job_name = "#PBS -N {}\n".format(ifunc);
                fid_i.write(job_name);
                
                err_file = "#PBS -e {}/{}.err\n".format(PBS_HOME, ifunc);
                fid_i.write(err_file);
                
                log_file = "#PBS -o {}/{}.log\n".format(PBS_HOME, ifunc);
                fid_i.write(log_file);

                move_workdir = "cd {}\n".format(ANALYZER_HOME);
                fid_i.write(move_workdir);
                
                run_job = "{0} {1}.py {2}\n".format(path_python, ifunc, file_out);
                fid_i.write(run_job);
                fid_i.close();
                
            # ---------------------------- Writing SHELL SCRIPT
            for ifunc in func_name:
                tmp = "{}/{}.sh".format(SH_HOME, ifunc);
                fid_i = open(tmp, 'w');
                
                def_shell = "#!/bin/bash\n";
                fid_i.write(def_shell);

                move_workdir = "cd {}\n".format(ANALYZER_HOME);
                fid_i.write(move_workdir);
                
                run_job = "{0} {1}.py {2}\n".format(path_python, ifunc, file_out);
                fid_i.write(run_job);
                fid_i.close();
                
                # change script permission to 755
                os.chmod(tmp, 0755);
            
            
            # ---------------------------- Insert job into a table ---------------------------------
            #print "HERE works"
            conn = sqlite3.connect(dbName);
            c = conn.cursor();
            
            funcNames = ",".join(func_name);
            stime     = datetime.now().strftime("%Y-%m-%d %H:%M:%S");
            query = """INSERT INTO gui_job (name, proj_id, anaz, status, output, stime, etime) \
                    VALUES ("{0}", {1}, "{2}", "{3}", "{4}", "{5}", "{6}")""".format(job_title, pkey, funcNames, 'SENT', OUTPUT_HOME, stime, 'N/A');
            #print query
            c.execute(query);
            conn.commit();
            #print "Job table updated!"
            
            #print "# ---------------------------- Insert parameter into a table ---------------------------------"
            #print "Paras: {}".format(Paras);
            #print "func_name: {}".format(func_name);
            #print "ParaInfo: {}".format(ParaInfo);
            #calculating the most recent Job primary key
            query = """SELECT id FROM gui_job WHERE proj_id = {0} AND stime = "{1}" """.format(pkey, stime);
            #print query
            c.execute(query);
            job_pkey = c.fetchone();
            para_pkeys = [];
            #print job_pkey[0];
            for i in range(len(func_name)):
                para_pkey = [];
                for j in range(len(Paras[i])):
                    query = """INSERT INTO gui_parameter (job_id, anaz, para, val, status) \
                            VALUES ({0}, "{1}", "{2}", "{3}", "{4}")""".format(job_pkey[0], func_name[i], ParaInfo[i][j], Paras[i][j], 'SENT');
                    #print query;
                    c.execute(query);
                    conn.commit();
                query = """SELECT id FROM gui_parameter WHERE job_id = {0} AND anaz = "{1}" """.format(job_pkey[0], func_name[i]);
                #print query
                c.execute(query);
                tmp_para_pkey = c.fetchall();
                for item in tmp_para_pkey:
                    para_pkey.append(item[0]);
                para_pkeys.append(para_pkey);
                #print "parameter table updated!"
                    
            #print "# ---------------------------- Display tables ---------------------------------"
            #print "Items in gui_job"
            query = "SELECT * from gui_job";
            c.execute(query);
            row = c.fetchall();
            #for item in row:
            #    print item;
                
            #print "Items in gui_parameter"
            query = "SELECT * from gui_parameter";
            c.execute(query);
            row = c.fetchall();
            #for item in row:
            #    print item;

            conn.close();


            #-------------------------------------------
            # Writring input arguments into a file
            #-------------------------------------------
            c ={
                'para_pkeys'        : para_pkeys,
                'job_pkey'          : job_pkey,
                'job_title'         : job_title,
                'pkey'              : pkey,
                'ptitle'            : ptitle,
                'dbName'            : dbName,
                'pbs'               : pbs,
                'date'              : date,
                'session_home'      : SESSION_HOME,
                'output_home'       : OUTPUT_HOME,
                'pbs_home'          : PBS_HOME,
                'media_home'        : MEDIA_HOME,
                'analyzer_home'     : ANALYZER_HOME,
                'base_path'         : bpath,
                'path_output'       : path_output,
                'path_python'       : path_python,
                'structure_file'    : stfile,
                'trajectory'        : trjFile,
                'Paras'             : Paras,
                'funcName'          : func_name,
                'paraInfo'          : ParaInfo,
            }
            
            #print "works {}".format(PBS_HOME)
            # save dictionary into a file
            #file_out = "{}para".format(OUTPUT_HOME)
            #print "Writing binary file {}".format(file_out);
            fid_out = open(file_out, 'wb');
            pickle.dump(c, fid_out);
            fid_out.close();
            
            # ---------------------------- Submit Jobs
            if (machine == 'cluster'):
                print "Run code at {}".format(machine);
                for ifunc in func_name:
                    cmd = "qsub {}/{}.pbs".format(PBS_HOME, ifunc);
                    os.system(cmd);

            elif (machine == 'server'):
                print "Run code at {}".format(machine);
                for ifunc in func_name:
                    cmd = "{}/{}.sh".format(SH_HOME, ifunc);
                    os.system(cmd);

            return HttpResponse(json.dumps(c));


def makeDownload(request):
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )
    
    request.session.set_expiry(SESSION_TIME_OUT);
    output_id 	= request.GET.get('id');
    dfield = request.GET.get('dformat');
    conn = sqlite3.connect(dbName);
    c = conn.cursor();
    query = """select {0} from gui_outputs where id = {1}""".format(dfield, output_id);
    #print query;
    c.execute(query);
    row = c.fetchone();
    #print "the path of [{0}] in output_id ({1}) is {2}".format(dfield, output_id, row[0]);
    path = row[0];
    
    # link the file
    wrapper = FileWrapper( open( path, "r" ) )
    content_type = mimetypes.guess_type( path )[0]

    response = HttpResponse(wrapper, content_type = content_type)
    response['Content-Length'] = os.path.getsize( path ) # not FileField instance
    response['Content-Disposition'] = 'attachment; filename={0}'.format(smart_str( os.path.basename( path ) ));

    return response

def mediaLink(request, file_name):
    #print "=== I am in mediaLink ==="
    path = '/{}'.format(file_name);
    #print path
    wrapper = FileWrapper( open( path, "r" ) )
    content_type = mimetypes.guess_type( path )[0]

    response = HttpResponse(wrapper, content_type = content_type)
    response['Content-Length'] = os.path.getsize( path ) # not FileField instance
    response['Content-Disposition'] = 'attachment; filename={0}'.format(smart_str( os.path.basename( path ) ));

    return response

def wysFileManager(request):
    #print request
    c ={
            'output'        : 'okay updated!'
	}
    return HttpResponse(json.dumps(c));

def resultView_DBmanager(request):
    print "* I am in resultView_DBmanager";
    # check out authority 
    if 'user_id' not in request.session:
        c = {
                    'errMsg'	    : 'Session has been expired!',
                }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request) )

    request.session.set_expiry(SESSION_TIME_OUT);
    user_id = request.session['user_id'];
    
    if request.is_ajax() and (request.method == 'POST'):
        cmd   = request.POST.get('cmd');
    else:
        cmd = 'http';           # connection with URL
        c = {
                    'errMsg'	    : 'Session has been expired!',
            }
        template = 'gui/login.html';
        return render_to_response(template, c, context_instance = RequestContext(request))
    
    if cmd == 'delete':
        print "OKAY I am in delete DB manager (Result View)"
	table 	= request.POST.get('table');
	tmpIDs	= request.POST.get('IDs');
	flag_del= request.POST.get('del');
	
	# parsing IDs
	listIDs = tmpIDs.split(',');
	IDs = [];
	for num in listIDs:
	    if '-' in num:
		print "FOUND conataining '-': {}".format(num);
		tmp = num.split('-');
		num1 = int(tmp[0]);
		num2 = int(tmp[1]);
		if num1 < num2:
		    tmp2 = range(num1,num2+1);
		else:
		    tmp2 = range(num2,num1+1);
		for i in tmp2:
		    IDs.append(i);
	    else:
		tmp2 = int(num);
		IDs.append(tmp2);
	IDs = set(IDs);
	IDs = list(IDs);
	
	# deleting gui_project
	if table == 'gui_project':
	    print "** DELETE gui_project **";
	    delProjects(IDs, dbName, flag_del);
	elif table == 'gui_outputs':
	    print "** DELETE gui_outputs **";
	    delOutputs(IDs, dbName, flag_del);
	    
        c = {
            'fUsers'    : IDs,
        }
    return HttpResponse(json.dumps(c));

