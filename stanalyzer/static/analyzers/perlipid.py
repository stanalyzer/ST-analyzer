#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# perlipid.py
# Date: June 20, 2013
#
# purpose: calculating perlipid area based on  Voronoi diagram and Delaunay triangulation
#          Voronoi diagram may be releatively more favorable for identifying area of individual lipids or amino acids
#           Delaunay triangulation may be relatively more favorable for calculating surface area consisting of multiple atoms
#
# Programmer: Jong Cheol Jeong (mejcjeong@gamil.com)
#
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# for local functions and classes
from datetime import datetime
import os
import subprocess
import stanalyzer


#---------------------------------------------------------------------
# Insert points to the maximum and minimum boundary
# of the system box
#---------------------------------------------------------------------
def fill2d3d(sysSize):
    for i in range(len(sysSize)):
        sysSize[i] = float(sysSize[i]);
    
    sidx = 0;
    sysMtx = [];
    
    sysMtx.append(sysSize);
    
    tmp = sysSize[:]; # copy list
    for i in range(len(sysSize)):
        tmp[i] = tmp[i] * -1;
    
    sysMtx.append(tmp);
    
    for i in range(len(sysSize)):
        tmp = sysSize[:]; # copy list
        tmp[i] = tmp[i] * -1;
        sysMtx.append(tmp);
    
    if len(sysSize) > 2:    
        for i in range(len(sysSize)):
            tmp = sysSize[:]; # copy list
            if i < (len(sysSize) -1):
                for j in frange(i,i+1, 1):
                    tmp[int(j)] = tmp[int(j)] * -1;
            else:
                tmp[i] = tmp[i] * -1;
                tmp[0] = tmp[0] * -1;
                
            sysMtx.append(tmp);
    return sysMtx;




#---------------------------------------------------------------------------------------------------
# calculating surface area for the given Vertices
# by using Voronoi diagram
# *** Parameters ***
# qhull_path  : location where qhull is installed (e.g. /usr/bin/qhull/)
#               qhull is available at http://www.qhull.org/
# CRD         : list type of coordinates CRD = [[x0,y0,z0], [x1,y1,z1],...,[xn,yn,zn],]
# sysSize     : system size. sysSize = [max_x, max_y];
# TMP_DIR    : temporary directory where qhull creates its output.
#          output files will be automatically deleted when the processing is done
#
# *** Return values ****
# areaDic['area']        # area calculated by voronoi diagram
# areaDic['system']      # area calculated by system box size
# areaDic['region']      # individual voronoi regions
# areaDic['coord']       # coordinates of input data
# areaDic['vertex']      # all vertices consisting of voronoi region
# areaDic['refidx']      # reference of vertex index for consisting each voronoi region
#
#---------------------------------------------------------------------------------------------------
def perlipidVro (TMP_DIR, qhull_path, CRD, sysSize):

    """------------------------------------------------------
    # write coordinates into a file
    # 1st line : dimension
    # 2nd line : number of vertex
    # rest lines: coordinates x, y, z
    ---------------------------------------------------------"""
    stime = datetime.now().strftime("%Y%m%d%H%M%S%f");
    rnum  = stanalyzer.rand_N_digits(5);
    #crdFile = "{0}/crd_{1}{2}.dat".format(TMP_DIR, stime, rnum);
    crdFile = "{0}/crd.dat".format(TMP_DIR);
    
    fcrd_id = open(crdFile, "w");
    dim = 2;                           # dimension - currently 2D only
    num_vtx = len(CRD) + 4;            # adding vertices at maximum & minimum corner
    
    sline = "{0}\n{1}\n".format(dim, num_vtx);
    fcrd_id.write(sline)

    line_cnt = 0;
    for vi in CRD:
        sline = "{0} {1}\n".format(vi[0], vi[1]);  # considering only X and Y axis
        fcrd_id.write(sline);

    # filling out four vertices at the minimum and maximum system box
    sysMtx = fill2d3d(sysSize);
    for bsys in sysMtx:
        lastLine = "";
        for baxi in bsys:
            tmp = abs(baxi) - 0.01;
            if baxi < 0:
                tmp = tmp * -1;
            lastLine = "{0} {1}".format(lastLine, tmp);
        lastLine = "{0}\n".format(lastLine);
        fcrd_id.write(lastLine);

    fcrd_id.close();


    """------------------------------------------------------
    # Run qhull to get Voronoi vertices and cell region
    ---------------------------------------------------------"""
    # qhull input file = coordinate file created above

    # qhull output file - all unique vertices for voronoi region
    stime = datetime.now().strftime("%Y%m%d%H%M%S%f");
    rnum  = stanalyzer.rand_N_digits(5);
    #vtxFile = "{0}/vtx_{1}{2}.dat".format(TMP_DIR, stime, rnum);
    vtxFile = "{0}/vtx.dat".format(TMP_DIR);
    
    # qhull output file - Vertices corresponding to each voronoi region
    stime = datetime.now().strftime("%Y%m%d%H%M%S%f");
    rnum  = stanalyzer.rand_N_digits(5);
    #rgnFile = "{0}/rgn_{1}{2}.dat".format(TMP_DIR, stime, rnum);
    rgnFile = "{0}/rgn.dat".format(TMP_DIR);
    
    qhull_path = stanalyzer.eval_path(qhull_path);
    qvoronoi="{0}qvoronoi".format(qhull_path);
    
    # writing shell script for running qhull
    shFile = "{0}/qhull.sh".format(TMP_DIR);
    qsh_id = open(shFile, "w");
    msg = "#!/bin/bash\n{0} p < {1} > {2}\n".format(qvoronoi, crdFile, vtxFile);
    msg = "{0}{1} FN < {2} > {3}\n".format(msg, qvoronoi, crdFile, rgnFile);
    qsh_id.write(msg);
    qsh_id.close();
    
    # chmode shell script
    qopt = ['chmod', '755', shFile];
    qhull = subprocess.call(qopt);
    
    # run shell script
    qopt= [shFile];
    qhull = subprocess.call(qopt);
    """
    # creating voronoi region
    qstring = "p TI {0} TO {1} ".format(crdFile, vtxFile);
    qopt = [qvoronoi, qstring];
    #qhull = subprocess.Popen(qopt, stdin=subprocess.PIPE, stdout=subprocess.PIPE);
    qhull = subprocess.call(qopt);
    
    # counting voronoi region
    qstring = "FN TI {0} TO  {1} ".format(crdFile, rgnFile);
    qopt = [qvoronoi, qstring];
    #qopt = [qvoronoi, "FN", "TI", crdFile, "TO", rgnFile];
    #qhull = subprocess.Popen(qopt, stdin=subprocess.PIPE, stdout=subprocess.PIPE);
    qhull = subprocess.call(qopt);
    """
    
    """------------------------------------------------------
    # Read vtxFile and rgnFile and save them in to list
    ---------------------------------------------------------"""
    Vtex = [];          # voronoi vertices
    Refi = [];          # voronoi region reference index

    # --- read vtxFile
    # 1st line: dimension
    # 2nd line: number of vertices
    # rest    : coordinates
    fid_v = open(vtxFile, 'r');
    line_cnt = 0;
    for readLine in fid_v:
        #print readLine
        line_cnt = line_cnt + 1;
        if (line_cnt == 1):
            dim = map(int, readLine.strip().split());
            dim = dim[0];
            
        if (line_cnt == 2):
            num_vtx = map(int, readLine.strip().split());
            num_vtx = num_vtx[0];
        
        if (line_cnt > 2):
            #print readLine.strip().split()
            tmp = map(float, readLine.strip().split());
            tmp = tmp[:dim];                                # extract coordinate based one the given dimension
            Vtex.append(tmp);
    fid_v.close();


    # --- read rgnFile
    # 1st line: number of voronoi regions
    # rest    : references [#of vertices at this region, list of vertex1,...n]
    fid_r = open(rgnFile, 'r');
    line_cnt = 0;
    for readLine in fid_r:
        line_cnt = line_cnt + 1;
        if (line_cnt == 1):
            num_ref = map(int, readLine.strip().split());
            
        if (line_cnt > 1):
            tmp = map(int, readLine.strip().split());
            Refi.append(tmp);
    fid_r.close();


    """------------------------------------------------------
    # Calculating the area of voronoi diagram
    ---------------------------------------------------------"""
    #print "- Total # of voronoi vertices: {}".format(len(Vtex));
    #print "- Total # of voronoi regions: {}".format(num_ref);
    
    #print "filtering vertices... "

    # List after filtering vertices located outside of system boundary
    newRefi = [];

    # valid Refi row index
    # vreg = [#of vertices corresponding to this region,
    #        index of vertices corresponding to the order of coordinates in vtxFile]
    for vreg in Refi:
        vregIdx = range(1, vreg[0]+1);
        tmpR = [];
        tmpR.append(0);     # number of ref. vertices
        for vidx in vregIdx:
            # if reference index < 0 : out of voronoi region
            if (vreg[vidx]) >= 0 :
                # vertices coordinates should be located inside the maximum and minimum system boundary
                if (abs(Vtex[vreg[vidx]][0]) <= float(sysSize[0])) and (abs(Vtex[vreg[vidx]][1]) <= float(sysSize[1])) :
                   tmpR.append(vreg[vidx]);
    
        # recounting the number of valid vertices for the current voronoi region
        tmpR[0] = len(tmpR)-1;
        # Collecting only valid region
        if (tmpR[0] > 0):
            newRefi.append(tmpR);
    #print "- The number of Voronoi region has been reduced: {0} -> {1}".format(len(Refi), len(newRefi));


    """------------------------------------------------------
    # Actual calculation of area
    ---------------------------------------------------------"""
    Varea = 0.0;            # total area calculated by Voronoi region
    Sarea = 0.0;            # total area calculated by system size
    tmpArea = [];           # area for each Voronoi region
    
    for vreg in newRefi:
        vregIdx = range(1, len(vreg));
        tmp_area = 0.0;
        for vidx in vregIdx:
            if vidx < vregIdx[len(vregIdx)-1]:
                x0 = Vtex[vreg[vidx]][0];
                y0 = Vtex[vreg[vidx]][1];
                x1 = Vtex[vreg[vidx+1]][0];
                y1 = Vtex[vreg[vidx+1]][1];
            else:
                x0 = Vtex[vreg[vidx]][0];
                y0 = Vtex[vreg[vidx]][1];
                x1 = Vtex[vreg[1]][0];
                y1 = Vtex[vreg[1]][1];
                
            tmp_area = tmp_area + 0.5 * (x0 * y1 - y0 * x1);
            tmpArea.append(tmp_area);
            
        Varea = Varea + abs(tmp_area);
    Sarea = 4 * sysSize[0] * sysSize[1];
 
    areaDic = {
            'area'      : Varea,        # area calculated by voronoi diagram
            'system'    : Sarea,        # area calculated by system box size
            'region'    : tmpArea,      # individual voronoi regions
            'coord'     : CRD,          # coordinates of input data
            'vertex'    : Vtex,         # all vertices consisting of voronoi region
            'refidx'    : newRefi       # reference of vertex index for consisting each voronoi region
            };
    
    # remove temporary files
    os.remove(crdFile);
    os.remove(vtxFile);
    os.remove(rgnFile);
    
    return areaDic;


#---------------------------------------------------------------------------------------------------
# calculating surface area for the given Vertices
# by using Delaunay triangulation
# *** Parameters ***
# qhull_path  : location where qhull is installed (can be downloaded from http://www.qhull.org/)
# CRD         : list type of coordinates CRD = [[x0,y0,z0], [x1,y1,z1],...,[xn,yn,zn],]
# sysSize     : system size. sysSize = [max_x, max_y];
# TMP_DIR    : temporary directory where qhull creates its output.
#          output files will be automatically deleted when the processing is done
#
# *** Return values ****
# areaDic['area']        # area calculated by delaunay triangulation
# areaDic['system']      # area calculated by system box size
# areaDic['region']      # individual delaunay triangulation regions
# areaDic['coord']       # coordinates of input data
# areaDic['vertex']      # all vertices consisting of delaunay triangulation which is same as areaDic['coord'];
# areaDic['refidx']      # reference of vertex index for consisting each delaunay triangulation
#
#---------------------------------------------------------------------------------------------------
def perlipidDT (TMP_DIR, qhull_path, CRD, sysSize):

    """------------------------------------------------------
    # write coordinates into a file
    # 1st line : dimension
    # 2nd line : number of vertex
    # rest lines: coordinates x, y, z
    ---------------------------------------------------------"""
    stime = datetime.now().strftime("%Y%m%d%H%M%S%f");
    rnum  = stanalyzer.rand_N_digits(5);
    #crdFile = "{0}/dt_crd_{1}{2}.dat".format(TMP_DIR, stime, rnum);
    crdFile = "{0}/dt_crd.dat".format(TMP_DIR);
    
    fcrd_id = open(crdFile, "w");
    dim = 2;                           # dimension - currently 2D only
    num_vtx = len(CRD);
    
    sline = "{0}\n{1}\n".format(dim, num_vtx);
    fcrd_id.write(sline)

    line_cnt = 0;
    for vi in CRD:
        sline = "{0} {1}\n".format(vi[0], vi[1]);  # considering only X and Y axis
        fcrd_id.write(sline);
        
    fcrd_id.close();


    """------------------------------------------------------
    # Run qhull to get the region of Delaunay triangulation 
    ---------------------------------------------------------"""
    # qhull input file = coordinate file created above

    # qhull output file - Vertices corresponding to each delaunay region
    stime = datetime.now().strftime("%Y%m%d%H%M%S%f");
    rnum  = stanalyzer.rand_N_digits(5);
    #rgnFile = "{0}/dt_rgn_{1}{2}.dat".format(TMP_DIR, stime, rnum);
    rgnFile = "{0}/dt_rgn.dat".format(TMP_DIR);
    
    qhull_path = stanalyzer.eval_path(qhull_path);
    qdelaunay = "{0}qdelaunay".format(qhull_path);


    # writing shell script for running qhull
    shFile = "{0}/qhull_dt.sh".format(TMP_DIR);
    qsh_id = open(shFile, "w");
    msg = "#!/bin/bash\n{0} i < {1} > {2}\n".format(qdelaunay, crdFile, rgnFile);
    qsh_id.write(msg);
    qsh_id.close();
    
    # chmode shell script
    qopt = ['chmod', '755', shFile];
    qhull = subprocess.call(qopt);
    
    # run shell script
    qopt= [shFile];
    qhull = subprocess.call(qopt);

    # counting delaunay triangulation region
    #qopt = [qdelaunay, "i", "TI", crdFile, "TO", rgnFile];
    #subprocess.Popen(qopt).communicate();

    """------------------------------------------------------
    # rgnFile and save them in to list
    ---------------------------------------------------------"""
    Vtex = CRD;         # vertices are now same as coordinates
    Refi = [];          # delaunay triangulation region reference index

    # --- read rgnFile
    # 1st line: number of voronoi regions
    # rest    : references [#of vertices at this region, list of vertex1,...n]
    fid_r = open(rgnFile, 'r');
    line_cnt = 0;
    for readLine in fid_r:
        line_cnt = line_cnt + 1;
        if (line_cnt == 1):
            num_ref = map(int, readLine.strip().split());
            
        if (line_cnt > 1):
            tmp = map(int, readLine.strip().split());
            Refi.append(tmp);
    fid_r.close();


    """------------------------------------------------------
    # Calculating the area of delaunay triangulation
    ---------------------------------------------------------"""
    #print "- Total # of voronoi vertices: {}".format(len(Vtex));
    #print "- Total # of voronoi regions: {}".format(num_ref);
    
    Darea = 0.0;            # total area calculated by delaunay triangulation
    Sarea = 0.0;            # total area calculated by system size
    tmpArea = [];           # area for each delaunay triangulation
    
    for vreg in Refi:
        tmp_area = 0.0;
        for vidx in range(len(vreg)):
            if vidx < (len(vreg)-1):
                x0 = Vtex[vreg[vidx]][0];
                y0 = Vtex[vreg[vidx]][1];
                x1 = Vtex[vreg[vidx+1]][0];
                y1 = Vtex[vreg[vidx+1]][1];
            else:
                x0 = Vtex[vreg[vidx]][0];
                y0 = Vtex[vreg[vidx]][1];
                x1 = Vtex[vreg[0]][0];
                y1 = Vtex[vreg[0]][1];
                
            tmp_area = tmp_area + 0.5 * (x0 * y1 - y0 * x1);
            tmpArea.append(tmp_area);
            
        Darea = Darea + abs(tmp_area);
    Sarea = 4 * sysSize[0] * sysSize[1];
 
    areaDic = {
            'area'      : Darea,        # area calcualted by delaunay triangulation
            'system'    : Sarea,        # area calculated by system box size
            'region'    : tmpArea,      # individual delaunay triangle regions
            'coord'     : CRD,          # coordinates of input data
            'vertex'    : Vtex,         # all vertices consisting of delaunay triangulation region (i.e. it is equal to CRD)
            'refidx'    : Refi          # reference of vertex index for consisting each triangle region  
            };
    
    # remove temporary files
    os.remove(crdFile);
    os.remove(rgnFile);
    
    return areaDic;
