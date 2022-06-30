import numpy as np
from numpy.linalg import eig
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
import math
plt.switch_backend('agg')

#alignment function - mostly identical to MDAnalysis.align._fit_to, but removes unnecessary
#translation step, still could incorporate bead weight
#see documentation at https://docs.mdanalysis.org/stable/_modules/MDAnalysis/analysis/align.html
def modFitTo(mobile_coordinates, ref_coordinates, mobile_atoms, mobile_com, weights=None):
    R, min_rmsd = align.rotation_matrix(mobile_coordinates, ref_coordinates, weights=weights)

    mobile_atoms.translate(-mobile_com)
    mobile_atoms.rotate(R)

    return mobile_atoms, min_rmsd

#full aligner function - aligns every monomer in a universe iterating by the numbers
#in monomer and time strides - needs to return an array containing all select positional data
def monAlign(universe, timeStride, monStride):
    #get centre of mass + position data for alignment reference
    #based on monomer halfway up the simulated section
    monHalf = math.trunc(len(universe.atoms.split("residue"))/2)
    refCOM = universe.atoms.split("residue")[monHalf].center_of_mass() #call method
    refPos = universe.atoms.split("residue")[monHalf].positions[:] - refCOM[:] #subtract COM
    
    #iterate with while loops so odd numbers don't break the code
    t = 0
    
    #prepare positions list (Coordinates Over Time)
    global COT
    COT = list()
    
    
    while t < len(u.trajectory):
        
        universe.trajectory[t]
        
        m = 0
        
        while m < len(u.atoms.split("residue")):
            
            #in each case monomer in question must be aligned
            mob = universe.atoms.split("residue")[m]
            mobCOM = mob.center_of_mass()
            mobPos = mob.positions[:] - mobCOM[:]
            
            mob, rmsd = modFitTo(mobPos, refPos, mob, mobCOM, weights = None)
            
            #LOOK UP VSTACK (might be more efficient than this)
            COT.append(mob.positions.tolist())
            
            #finish monomer operations
            m += monStride
        #finish frame operations
        t += timeStride
    
    #return Coordinates over time
    COTarray = np.asarray(COT)
    return COTarray

#monalign returns data in format: (monomer, bead, coordinate)
#bead and coordinate may need to be combined
def combine_dims(a, start=0, count=2):
    s = a.shape
    return np.reshape(a, s[:start] + (-1,) + s[start+count:])

#x and y must each be vectors that can be iterated through
def CoV(x, y):
    ex = np.mean(x)
    ey = np.mean(y)
    
    sumCov = 0
    
    for xj, yj in zip(x, y):
        ix = xj - ex
        iy = yj - ey
        
        sumCov += ix*iy
        
    CoV = sumCov / len(x)
    return CoV

#makes a covariance matrix based on the coordinates - has to account for structure of the 
#arrays (monomer, bead, coordinate)
def monomerCoVarMatrix(coords):
    
    mxSize = len(coords[0,:])*3
    CoVMatrix = np.zeros([mxSize, mxSize])
    
    mxx = 0
    for beadx in range(coords.shape[1]):
        for xyzx in range(coords.shape[2]):
            
            mxy = 0
            for beady in range(coords.shape[1]):
                for xyzy in range(coords.shape[2]):
                    
                    CoVMatrix[mxx,mxy] = CoV(coords[:,beadx,xyzx], coords[:,beady,xyzy])
                    mxy += 1
            mxx += 1
    return CoVMatrix


def monomerPCA(coords):
    cvMatrix = monomerCoVarMatrix(coords)
    eigvalues, eigvectors = eig(cvMatrix)
    PCs = np.dot(combine_dims(coords), eigvectors)
    
    return PCs, eigvalues

#polymer and run names for file access
polys = list(["5o4u", "pilA4", "TTC1836"])
runs = list(["3", "2", "7", "4", "5"])

#stride values for reducing resolution to save time, iterate by appropriate stride value
tStd = 100 #stride to count frames
mStd = 8   #stride to count monomers

#for every polymer type, and every run, load in universe
for poly in polys:
    for run in runs:
        u = mda.Universe("{}_cg.psf".format(poly), "{}_Run_{}.dcd".format(poly,run))
        
        coords = monAlign(u, tStd, mStd)
        PCA, eigenvalues = monomerPCA(coords)
