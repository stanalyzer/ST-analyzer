def perlipidVro (funcName, dic):
    

psf = "/Users/jcjeong/project/trajectory/step5_assembly.psf";
pdb = "/Users/jcjeong/project/trajectory/step7_691.pdb";
dcd = "/Users/jcjeong/project/trajectory/step7_633.dcd";

MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False

u = Universe(psf,dcd);
qry = "segid MEMB and resname DOPC and name C2* and not (name C2 or name C21)";

selAtoms = u.selectAtoms(qry);

#---------------------------------------------------
# parameters
#---------------------------------------------------
qhull_path = "/usr/local/bin/"

sysSize = [0.5, 0.5];

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



# get coordinates
#CRD = selAtoms.coordinates();

CRD = [[-0.5, -0.3685, 0.1], 
        [-0.453, 0.1789, 0.1],
        [0.453, -0.1789, 0.1],
        [-0.453, 0.2789, 0.1],
        [-0.453, -0.3789, 0.1],
        [-0.453, -0.4789, 0.1],
        [-0.453, 0.3789, 0.1],
        [0.453, -0.2789, 0.1],
        [-0.353, 0.1789, 0.1],
        [0.253, -0.0789, 0.1],
        [-0.153, 0.2789, 0.1],
        [-0.053, -0.3789, 0.1],
        [-0.153, -0.4789, 0.1],
        [-0.253, 0.3789, 0.1],
        [0.333, -0.2789, 0.1],
        [-0.453, 0.1789, 0.1],
        [0.353, -0.0789, 0.1],
        [0.253, 0.1789, 0.1],
        [0.153, -0.2789, 0.1],
        [-0.4465, 0.0297, 0.1],
        [-0.1166, -0.4332, 2.0],
        [-0.1165, 0.0194, 0.1], 
        [0.0328, -0.281, 0.1], 
        [0.1711, -0.4923, 1.0],
        [0.1793, 0.4347, 0.1],
        [0.2556, -0.0413, 0.1],
        [0.331, -0.4654, 0.1], 
        [0.491, -0.4954, 0.1],
        [0.491, 0.4954, 0.1],
        [-0.491, -0.4954, 0.1],
        [0.491, 0.4954, 0.1]];

#---------------------------------------------------
# write coordinates into a file
# 1st line : dimension
# 2nd line : number of vertex
# rest lines: coordinates x, y, z
#---------------------------------------------------
out_fname = "xyz.dat";
f = open("xyz.dat", "w");
dim = 2;
num_vtx = len(CRD) + 4;
sline = "{0}\n{1}\n".format(dim, num_vtx);
f.write(sline)

line_cnt = 0;
for vi in CRD:
    sline = "{0} {1}\n".format(vi[0], vi[1]);
    f.write(sline);


sysMtx = fill2d3d(sysSize);
for bsys in sysMtx:
    lastLine = "";
    for baxi in bsys:
        tmp = abs(baxi) - 0.01;
        if baxi < 0:
            tmp = tmp * -1;
        lastLine = "{0} {1}".format(lastLine, tmp);
    lastLine = "{0}\n".format(lastLine);
    f.write(lastLine);

f.close();


#---------------------------------------------------
# Run qhull to get Voronoi vertices and cell region
#---------------------------------------------------
qh_in = out_fname;
qhv_out = "v.dat";
qhr_out = "r.dat";

qvoronoi="{0}qvoronoi".format(qhull_path);

# calculating voronoi vertices
qopt = [qvoronoi, "p", "TI", qh_in, "TO", qhv_out];
subprocess.Popen(qopt).communicate();

# counting voronoi region
qopt = [qvoronoi, "FN", "TI", qh_in, "TO", qhr_out];
subprocess.Popen(qopt).communicate();

#---------------------------------------------------
# Read v.dat and r.dat and save them in to list
#---------------------------------------------------
Vtex = [];          # voronoi vertices
Refi = [];          # voronoi region reference index

# read file v.dat
# 1st line: dimension
# 2nd line: number of vertices
# rest    : coordinates
fid_v = open(qhv_out, 'r');
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

# read file r.dat
# 1st line: number of references
# rest    : references [#of vertices vertex1,...n]

fid_r = open(qhr_out, 'r');
line_cnt = 0;
for readLine in fid_r:
    line_cnt = line_cnt + 1;
    if (line_cnt == 1):
        num_ref = map(int, readLine.strip().split());
        
    if (line_cnt > 1):
        tmp = map(int, readLine.strip().split());
        Refi.append(tmp);

fid_r.close();


#---------------------------------------------------
# Calculating area 
#---------------------------------------------------
print "total # of voronoi vertices: {}".format(len(Vtex));

print "filtering vertices... "

newRefi = [];

# valid Refi row index
for vreg in Refi:
    print vreg
    vregIdx = range(1, vreg[0]+1);
    tmpR = [];
    tmpR.append(0);     # number of ref. vertices
    print vregIdx
    for vidx in vregIdx:
        # reference index is negative then out of voronoi region
        if (vreg[vidx]) >= 0 :
            """
            if (abs(Vtex[vreg[vidx]][0]) > float(sysSize[0])):
                Vtex[vreg[vidx]][0] = sysSize[0];
                
            if (abs(Vtex[vreg[vidx]][1]) > float(sysSize[1])):
                Vtex[vreg[vidx]][1] = sysSize[1];
                
            tmpR.append(vreg[vidx]);
            """
            if (abs(Vtex[vreg[vidx]][0]) <= float(sysSize[0])) and (abs(Vtex[vreg[vidx]][1]) <= float(sysSize[1])) :
               tmpR.append(vreg[vidx]);

    tmpR[0] = len(tmpR)-1;
    #print tmpR
    if (tmpR[0] > 0):
        newRefi.append(tmpR);

print "Refi has been reduced: {0} -> {1}".format(len(Refi), len(newRefi));

print "Calculating area... "
Tarea = 0.0;
Tarea2 = 0.0;

out_fname = "output.dat";
f = open(out_fname, "w");

for vreg in newRefi:
    print vreg
    #vregIdx = range(1, vreg[0]+1);
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
        
        print "x0={}, y0={}\tx1={}, y1={}".format(x0, y0, x1, y1);
        tmp_area = tmp_area + 0.5 * (x0 * y1 - y0 * x1);
    tmpLine = "{0}\n".format(abs(tmp_area));
    print "tmp Area = {}".format(tmp_area);
    f.write(tmpLine);
        
    Tarea = Tarea + abs(tmp_area);

sarea = (2* sysSize[0]) * (2 * sysSize[1]);
print "system area = {0}, lipid area = {1}".format(sarea, abs(Tarea));
f.close();

