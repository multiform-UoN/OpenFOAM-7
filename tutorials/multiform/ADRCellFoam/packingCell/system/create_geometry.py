from math import *
from numpy import *
import sys
import subprocess


print("Create packing geometry\n----------------")
print("Usage: create_geometry.py packing_type porosity ncells resolution")
print("Usage: last argument can always be a boolean to enable voxelized mesh")

def print_preamble(f):
    with open("geometry1","r") as f2:
        for line in f2:
            f.write(line)
    return


def print_postamble(f):
    f.write("\n}\n};\n")
    return


def unitcell(p):
    if p=="bcc":
        x = [0.5,0.]#,0.,0.,1.]#,1.,1.,1.,1.]
        y = [0.5,0.]#,0.,1.,0.]#,1.,1.,0.,0.]
        z = [0.5,0.]#,1.,0.,0.]#,1.,0.,1.,0.]
    elif p=="sc":
        x = [0.5]#,0.,0.,1.]#,1.,1.,1.,1.]
        y = [0.5]#,0.,1.,0.]#,1.,1.,0.,0.]
        z = [0.5]#,1.,0.,0.]#,1.,0.,1.,0.]
    elif p=="fcc":
        x = [0.5,0.,0.5,0.]#,0.,0.,1.]#,1.,1.,1.,1.]
        y = [0.5,0.5,0.,0.]#,0.,1.,0.]#,1.,1.,0.,0.]
        z = [0.,0.5,0.5,0.]#,1.,0.,0.]#,1.,0.,1.,0.]
    return x,y,z


def write_domainsize(ls,ly,lz):

    return


def write_unitcell(f,ptype,x0,y0,z0,l):
    a,b,c=unitcell(ptype)
    j=0
    for x,y,z in zip(a,b,c):
        xx=x*l+x0
        yy=y*l+y0
        zz=z*l+z0
        f.write("        grain%d\n" % (n+j))
        f.write("        {\n")
        f.write("                surface circle;\n")
        f.write("                scale   (0.5 0.5 0.5);\n")
        f.write("                transform\n")
        f.write("                {\n")
        f.write("                        coordinateSystem\n")
        f.write("                {\n")
        f.write("                        type    cartesian;\n")
        f.write("                        origin  (%f %f %f);\n" % (xx,yy,zz))
        f.write("                coordinateRotation\n")
        f.write("                {\n")
        f.write("                type    axesRotation;\n")
        f.write("                e1      (1 0 0);\n")
        f.write("                e3      (0 0 1);\n")
        f.write("                }\n")
        f.write("                }\n")
        f.write("                }\n")
        f.write("        }\n")
        j+=1


def unit_length(t,p):
    if t=="bcc":
        dens=1-p
        d1sphere=(dens/2.*6./math.pi)**(1./3.)
        return 1/d1sphere
    elif t=="sc":
        dens=1-p
        d1sphere=(dens*6./math.pi)**(1./3.)
        return 1/d1sphere
    elif t=="fcc":
        dens=1-p
        d1sphere=(dens/4.*6./math.pi)**(1./3.)
        return 1/d1sphere

def find_location(t,l):
    return 0.5,0.01,0


print(sys.argv)
ptype     = sys.argv[1] # packing type bcc, scc, fcc
porosity  = float(sys.argv[2])
ncells    = float(sys.argv[3])
res    = float(sys.argv[4])
try:
    ncelly = float(sys.argv[5])
    ncellz = float(sys.argv[6])
    ncellx = ncells
except:
    ncelly=ncells
    ncellx=ncells
    ncellz=ncells
snap=True

print(ptype,porosity,ncells,ncellx,ncelly,ncellz,res,snap)

if (sys.argv[-1] in ("T", "t", "true", "TRUE", "True")):
    snap=False

l = unit_length(ptype,porosity)
with open("geometry","w") as f:
    print_preamble(f)
    n=0
    for i in arange(ceil(ncellx)+1):
        for j in arange(ceil(ncelly)+1):
            for k in arange(ceil(ncellz)+1):
                write_unitcell(f,ptype,l*(i),l*(j),l*(k),l)
                n+=10
    print_postamble(f)

(xin,yin,zin)=find_location(ptype,l)
x1=0
y1=0
z1=0
x2=ncellx
y2=ncelly
z2=ncellz
with open("domainsize","w") as f:
    f.write("scalegrid %f;\n" % (l))
    f.write("x1 %f;\n" % (x1))
    f.write("y1 %f;\n" % (y1))
    f.write("z1 %f;\n" % (z1))
    f.write("xin %f;\n" % (xin))
    f.write("yin %f;\n" % (yin))
    f.write("zin %f;\n" % (zin))
    f.write("x2 %f;\n" % (x2))
    f.write("y2 %f;\n" % (y2))
    f.write("z2 %f;\n" % (z2))
    f.write("xgrid %d;\n" % (int(x2*res)))
    f.write("ygrid %d;\n" % (int(y2*res)))
    f.write("zgrid %d;\n" % (int(z2*res)))
    if snap:
        f.write("snap true;\n")
    else:
        f.write("snap false;\n")
    f.write("#inputMode merge")

print("DONE!")
