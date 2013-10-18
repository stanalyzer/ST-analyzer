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
from bisect import bisect_left


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

def crdArea(Region, Vtex, maxX, maxY):
    # Region: index of vertices building a region
    # Vtex  : vertices of regions
    # filetering region if coordinates exceeds maxY and maxX
    newRegion = [];
    for rgn in Region:
        x = Vtex[rgn][0];
        y = Vtex[rgn][1];
        if (x <= maxX) and (x >= -maxX) and (y <= maxY) and (y >= -maxY):
            newRegion.append(rgn);
    Region = newRegion;
    area = 0.0;    
    if len(newRegion) > 2:
        for idx in range(len(Region)):
            if idx >= (len(newRegion) -1):
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

def voroArea(topAtoms, qry, maxX, maxY):
    topCRDs = topAtoms.coordinates();
    
    #convert 3D to 2D: extractig X and Y axis
    topCRDs = selCoord(topCRDs, [0, 1]);
    topPoints = getNumPoints(topAtoms, qry);
    
    uq_resNames = topPoints[0];                    # Residues names in the given system
    uq_numAtoms = topPoints[1];                    # the number of selected atoms corresponding to residue names
    
    # calculating voronoi tessellation
    voro = VoronoiTess(topCRDs);
    voro_vtex = voro.vertices;
    voro_regn = voro.regions;
    
    # calculating area per lipids
    resIDs     = topAtoms.resids();               # obtaining residue ids
    resNames   = topAtoms.resnames();
    
    # initializing array 
    areaLipids = [0] * (len(resIDs)+1);           # each array will be used for individual lipid area and last one will be the average area
    
    res_cnt   = 0;                               # residue count
    num_atoms = 0;                               # number of atoms corresponding to the current residue name
    crd_idx   = 0;                               # coordinate index
    for res_name in resNames:
        # find corresponding residue
        num_atoms = uq_numAtoms[uq_resNames.index(res_name)];
        #calculating the area corresponding to the point
        atom_cnt = 0;
        myarea    = 0.0;
        for dummy in range(num_atoms):
            print crd_idx;
            my_region = getPos(voro_regn[crd_idx]);
            liparea = crdArea(my_region, voro_vtex, maxX, maxY);
            print liparea
            crd_idx = crd_idx + 1;
            myarea = myarea + liparea;
        areaLipids[res_cnt] = myarea;
        res_cnt = res_cnt + 1;
    areaLipids[len(areaLipids)-1] = sum(areaLipids) / res_cnt;
    print crd_idx
    print len(topCRDs)
    return areaLipids

def getCRDsWithResid(AtomGroup):
    CRDs = [];
    ATOMs = [];
    resIDs = AtomGroup.resids();
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
    resInfo = [resIDs, ATOMs, CRDs];	
    return resInfo

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

# convert degree to radian
def toRadian (fromDeg):
    rad = float(fromDeg) * math.pi / 180.0;
    return rad;

# convert radian to degree
def toDegree(fromRad):
    deg = float(fromRad) * 180.0 / math.pi;
    return deg;

#psf = "/home2/jcjeong/project/charmm/dol_dopc_1/step5_assembly.xplor_ext.psf";
psf = "/imscratch/yfqi/forJC/step1_pdbreader.psf";
#psf = "/home2/jcjeong/project/charmm/emilia2/yeast_pc_popi13_1/step5_assembly.psf";

#pdb = "/home2/jcjeong/project/charmm/dol_dopc_1/step5_assembly.pdb";
pdb = "/home2/jcjeong/project/charmm/emilia2/yeast_pc_popi13_1/step5_assembly.pdb";

#dcd = "/home2/jcjeong/project/charmm/dol_dopc_1/step7_100.dcd";
dcd = "/home2/jcjeong/project/charmm/emilia2/yeast_pc_popi13_1/step7_101.dcd"
dcd = "/imscratch/yfqi/forJC/step5_production.now.dcd";

MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False

u = Universe(psf,dcd);

qry = "segid MEMB and ((resname CHL1 and name O3) \
      or (resname DOPC and (name C2 or name C21 or name C31)) \
      or (resname POPC and (name C2 or name C21 or name C31)) \
      or (resname POPI and (name C2 or name C21 or name C31)) \
      or (resname POPI13 and (name C2 or name C21 or name C31)))";

selQry = "segid MEMB and resname CHL1 and (name C3 or name C17)";

#creating BIN 0~90 degree
BIN = []; DIST = [];
for ibin in stanalyzer.frange(0, 90, 0.1):
    BIN.append(ibin);
    DIST.append(0.0);


L = MDAnalysis.analysis.leaflet.LeafletFinder(u, selQry, cutoff=15.0, pbc=True);
topAtoms = L.group(0);
AtRes = getCRDsWithResid(topAtoms);
ResIDs  = AtRes[0]; # residue IDs [resid1, resid2, ..., residN]
AtmName = AtRes[1]; # [[name11, name12], [name21, name22],....,[nameN1, nameN2]] - assume each residue has 2 atoms
AtmCrds = AtRes[2] # [[[x11,y11,z11][x12,y12,z12]],[[x21,y21,z21][x22,y22,z22]] ...[[xN1,yN1,zN1][xN2,yN2,zN2]]] 
for ridx in range(len(ResIDs)):
    vtx1 = AtmCrds[ridx][0];	# atom C3
    vtx2 = AtmCrds[ridx][1];	# atom C17
    v0 = [0, 0, 1];		# bilayer normal
    v1 = [];
    for axi in range(len(vtx1)):
        tmp = vtx1[axi] - vtx2[axi];
        v1.append(tmp);
    T = toDegree(math.acos(getCosT(v0, v1)));
    pos = bisect_left(BIN, T);
    DIST[pos] = DIST[pos] + 1;

# box length
maxX = 50.0;
maxY = 50.0;

# coordinates
maxX = 0.5 * maxX;
maxY = 0.5 * maxY;
flg_top = 0;
flg_btm = 0;

# for bilayer membrane
topQry = "{0} and (prop z > 0.0)".format(qry);
topAtoms = u.selectAtoms(topQry);
if len(topAtoms) > 0:
    topArea = voroArea(topAtoms, qry, maxX, maxY);
    top_resIDs     = topAtoms.resids();               # obtaining residue ids
    flg_top = 1;

btmQry = "{0} and (prop z < 0.0)".format(qry);
btmAtoms = u.selectAtoms(btmQry);
if len(btmAtoms) > 0:
    btmArea = voroArea(btmAtoms, qry, maxX, maxY);
    btm_resIDs     = topAtoms.resids();               # obtaining residue ids
    btm_resNames   = topAtoms.resnames();
    flg_btm = 1;


# at the final stage calculating area per lipids
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
