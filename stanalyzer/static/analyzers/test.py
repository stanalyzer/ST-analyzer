import sys
from MDAnalysis import *
import datetime
from collections import defaultdict	# for initializing dictionary
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
from collections import defaultdict # for initializing dictionary

import stanalyzer

dnst_min = 2;
dnst_max = 18;
dnst_bin = 1;

BIN = [];
for ibin in frange(dnst_min, dnst_max, dnst_bin):
    BIN.append(ibin);

psf = "/home2/jcjeong/project/charmm/dol_dopc_1/step5_assembly.xplor_ext.psf";
pdb = "/home2/jcjeong/project/charmm/dol_dopc_1/step5_assembly.pdb";
dcd = "/home2/jcjeong/project/charmm/dol_dopc_1/step7_100.dcd";
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False
u = Universe(psf,dcd);
qry = "segid MEMB and resname DOPC and name C2* and not (name C2 or name C21)";

selAtoms = u.selectAtoms(qry);

Names = selAtoms.names();
setNames = set(Names);
listNames = list(setNames);

sortlistNames = stanalyzer.sort_str_num(listNames, "ASCE")
