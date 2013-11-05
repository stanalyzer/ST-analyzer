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
from pyhull.voronoi import VoronoiTess
import numpy


def getNumPoints(selAtoms, qry):
    numPoints = [];
    ResNames    = selAtoms.resnames();         # list of all residue names
    uqResNames  = list(set(ResNames));       # list of unique resiude names
    Names       = selAtoms.names();
    # parsing query to get information about the number of points per lipid
    uqResNames_num = [];
    byResname = qry.split('resname');
    for res in uqResNames:
        ptrn = "\\b{0}\\b".format(res);
        objptrn = re.compile(ptrn);
        for substr in byResname:
            tmp = substr.strip();
            #if res in tmp:
            mitem = objptrn.findall(tmp);
            if len(mitem) > 0:
                #print "{0} has {1} atoms".format(res, len(re.findall(r'name', tmp)));
                uqResNames_num.append(len(re.findall(r'name', tmp)));
    numPoints.append(uqResNames);
    numPoints.append(uqResNames_num);
    return numPoints

# this function extract coordinates based on the given index
# this function is used for creating 2D from 3D
def selCoord(CRDs, Idx):
    newCRDs =[];
    for cs in CRDs:
        tmp = [];
        for i in Idx:
            tmp.append(cs[i]);
        newCRDs.append(tmp);
    return newCRDs

# omitting negative index
def getPos(Rgn):
    newRgn = [];
    for i in Rgn:
        if i >= 0:
            newRgn.append(i);
    return newRgn

def crdArea(Region, Vtex, max_x, max_y, min_x, min_y):
    # Region: index of vertices building a region
    # Vtex  : vertices of regions
    # filetering region if coordinates exceeds maxY and maxX
    newRegion = [];
    for rgn in Region:
        x = Vtex[rgn][0];
        y = Vtex[rgn][1];
        if (x <= max_x) and (x >= min_x) and (y <= max_y) and (y >= min_y):
            newRegion.append(rgn);
    Region = newRegion;
    area = 0.0;    
    if len(Region) > 2:
        for idx in range(len(Region)):
            if idx >= (len(Region) -1):
                x1 = Vtex[Region[idx]][0];
                y1 = Vtex[Region[idx]][1];
                x2 = Vtex[Region[0]][0];
                y2 = Vtex[Region[0]][1];
            else:
                x1 = Vtex[Region[idx]][0];
                y1 = Vtex[Region[idx]][1];
                x2 = Vtex[Region[idx+1]][0];
                y2 = Vtex[Region[idx+1]][1];
                    
            area = area + 0.5 * (x1*y2 - x2*y1);
    return abs(area)

def mkimagebox(CRDs, size_x, size_y):
    # 5 |    3      | 7
    # 1 | primary 0 | 2
    # 6 |     4     | 8
    box0 = CRDs; box1 = []; box2 = [];
    box3 = [];   box4 = []; box5 = [];
    box6 = [];   box7 = []; box8 = [];

    for crd in CRDs:
        # box1
        x = crd[0] - size_x;
        y = crd[1];
        tmp = [x, y];
        box1.append(tmp);
        
        # box2
        x = crd[0] + size_x;
        y = crd[1];
        tmp = [x, y];
        box2.append(tmp);

        # box3
        x = crd[0];
        y = crd[1] + size_y;
        tmp = [x, y];
        box3.append(tmp);

        # box4
        x = crd[0];
        y = crd[1] - size_y;
        tmp = [x, y];
        box4.append(tmp);

        # box5
        x = crd[0] - size_x;
        y = crd[1] + size_y;
        tmp = [x, y];
        box5.append(tmp);

        # box6
        x = crd[0] - size_x;
        y = crd[1] - size_y;
        tmp = [x, y];
        box6.append(tmp);

        # box7
        x = crd[0] + size_x;
        y = crd[1] + size_y;
        tmp = [x, y];
        box7.append(tmp);

        # box8
        x = crd[0] + size_x;
        y = crd[1] - size_y;
        tmp = [x, y];
        box8.append(tmp);
        
    sysCRDs = box0 + box1 + box2 + box3 + box4 + box5 + box6 + box7 + box8;
    return sysCRDs

def voroArea(Atoms, qry, size_x, size_y):
    orgCRDs = Atoms.coordinates();
    
    #convert 3D to 2D: extractig X and Y axis
    orgCRDs = selCoord(orgCRDs, [0, 1]);
    orgPoints = getNumPoints(Atoms, qry);
    
    uq_resNames = orgPoints[0];                    # Residues names in the given system
    uq_numAtoms = orgPoints[1];                    # the number of selected atoms corresponding to residue names
    
    #-------------------------------
    # calculating image boxes
    # 5 |    3      | 7
    # 1 | primary 0 | 2
    # 6 |     4     | 8
    #-------------------------------
    sysCRDs = mkimagebox(orgCRDs, size_x, size_y);
    
    # calculating voronoi tessellation
    voro = VoronoiTess(sysCRDs);
    voro_vtex = voro.vertices;
    voro_regn = voro.regions;
    
    # calculating area per lipids
    resIDs     = Atoms.resids();               # obtaining residue ids
    resNames   = Atoms.resnames();
    
    # initializing array 
    areaLipids = [0.0] * (len(resIDs)+1);           # each array will be used for individual lipid area and last one will be the average area
    
    res_cnt   = 0;                               # residue count
    num_atoms = 0;                               # number of atoms corresponding to the current residue name
    crd_idx   = 0;                               # coordinate index
    M = numpy.array(sysCRDs);

    max_x = M.max(axis=0)[0] / 2.0;
    max_y = M.max(axis=0)[1] / 2.0;
    min_x = M.min(axis=0)[0] / 2.0;
    min_y = M.min(axis=0)[1] / 2.0;
    
    for res_name in resNames:
        # find corresponding residue
        num_atoms = uq_numAtoms[uq_resNames.index(res_name)];
        #print resNames;
        #calculating the area corresponding to the point
        atom_cnt = 0;
        myarea    = 0.0;
        for dummy in range(num_atoms):
            #print "max_x={}, max_y={}, min_x={}, min_y={}".format(max_x, max_y, min_x, min_y);
            #print "crdIdx:{0}, lenRegio:{1}, len(CRDs):{2}, residue:{3}, num_atoms:{4}, num_residue:{5}".format(crd_idx, len(voro_regn), len(orgCRDs), res_name, num_atoms, len(resNames));
            if crd_idx >= len(orgCRDs):
                break;
            else:
                my_region = getPos(voro_regn[crd_idx]);
                liparea = crdArea(my_region, voro_vtex, max_x, max_y, min_x, min_y);
                #print "crd_idx={}, res_name ={}, res_ID={}, X={}, Y={}".format(crd_idx, res_name, resIDs[res_cnt], sysCRDs[crd_idx][0], sysCRDs[crd_idx][1]);
                crd_idx = crd_idx + 1;
                myarea = myarea + liparea;
        areaLipids[res_cnt] = myarea;
        res_cnt = res_cnt + 1;
    areaLipids[len(areaLipids)-1] = sum(areaLipids) / res_cnt;
    #print "len(areaLipids): {}, len(orgCRDs): {}".format(len(areaLipids), len(orgCRDs));
    return areaLipids

def stateLipidArea(topAtoms, topArea):
    resInfo = [];
    top_resNames   = topAtoms.resnames();
    uq_resNames    = list(set(top_resNames));
    aveResArea =[];
    for res in uq_resNames:
        Idx = stanalyzer.findIndex(top_resNames, res);
        sp = Idx[0];
        ep = Idx[len(Idx)-1] + 1;
        ave = numpy.mean(topArea[sp:ep]);
        std = numpy.std(topArea[sp:ep]);
        stderr = std / math.sqrt(len(Idx));
        tmp = [ave, std, stderr];
        aveResArea.append(tmp);
    resInfo = [uq_resNames, aveResArea]
    return resInfo
