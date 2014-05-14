#!/opt/local/bin/python2.7

import os
#import numpy
import math
from numpy import *
from sys import *
import sys
import getopt

### --- Arguments --- ###
program = "XYZ-to-ForceField.py"
ifile = ''
ofile = ''

### Read command line args
try:
	myopts, args = getopt.getopt(sys.argv[1:],"i:o:h")
except getopt.GetoptError:
	print program + " -i <inputfile.map> -o <outputfile.map>"
	sys.exit(2)
###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
	if o == '-i':
		ifile = a
	elif o == '-o':
		ofile = a
	elif o == '-h':
		print program + " -i <inputfile.map> -o <outputfile.map>"
		sys.exit(0)
	else:
		print("Usage: %s -i inputfile.map" % sys.argv[0])
		sys.exit(0)
 
class Atom(object):
	def __init__(self):
		self.e = 0
		self.x = 0
		self.y = 0
		self.z = 0
		self.neighbors = []
		self.neighborsdist = []
		self.nearest = 0
		self.hybridization = ''
		self.charge = 0
		self.int = 0
		self.string = ''
		self.atomtype = ''


### --- Functions to get and give element numbers and names --- ###
elementList = ["h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p","s","cl","ar","k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te","i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb","lu","hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn","fr","ra","ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr","rf","db","sg","bh","hs","mt","ds","rg","cn","uut","fl","uup","lv","uus","uuo"]
elementNames = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"]
elementLarge =  ["na","mg","al","si","p","s","cl","ar","k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te","i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb","lu","hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn","fr","ra","ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr","rf","db","sg","bh","hs","mt","ds","rg","cn","uut","fl","uup","lv","uus","uuo"]
def getElementNum(at1):
	element = elementList.index(at1.lower())
	return element+1
def getElementName(at1):
	element = elementNames[at1-1]
	return element

### --- Distance function --- ###
def getdistance(at1, at2, lol):
	if isinstance(lol[at1], Atom):
		atom1 = array([lol[at1].x, lol[at1].y, lol[at1].z])
		atom2 = array([lol[at2].x, lol[at2].y, lol[at2].z])
	else:
		atom1 = array([float(lol[at1][1]), float(lol[at1][2]), float(lol[at1][3])])
		atom2 = array([float(lol[at2][1]), float(lol[at2][2]), float(lol[at2][3])])
	vector1 = atom2-atom1
	dist = linalg.norm(vector1)
	#print dist
	return dist

### --- Angle function --- ###
def getangle(at1, at2, at3, lol):
	if isinstance(lol[at1], Atom):
		# put positions in array
		atom1 = array([lol[at1].x, lol[at1].y, lol[at1].z])
		atom2 = array([lol[at2].x, lol[at2].y, lol[at2].z])
		atom3 = array([lol[at3].x, lol[at3].y, lol[at3].z])
	else:
		atom1 = array([float(lol[at1][1]), float(lol[at1][2]), float(lol[at1][3])])
		atom2 = array([float(lol[at2][1]), float(lol[at2][2]), float(lol[at2][3])])
		atom3 = array([float(lol[at3][1]), float(lol[at3][2]), float(lol[at3][3])])
	# making appropriate vectors and normals 
	vector1 = atom1-atom2
	vector2 = atom3-atom2
	angle = arccos(dot(vector1,vector2)/(linalg.norm(vector1)*linalg.norm(vector2)))
	#print degrees(angle)
	return degrees(angle)

### --- Dihedral angle function --- ###
def getdihedral(at1, at2, at3, at4, lol):
	if isinstance(lol[at1], Atom):
		# put positions in array
		atom1 = array([lol[at1].x, lol[at1].y, lol[at1].z])
		atom2 = array([lol[at2].x, lol[at2].y, lol[at2].z])
		atom3 = array([lol[at3].x, lol[at3].y, lol[at3].z])
		atom4 = array([lol[at4].x, lol[at4].y, lol[at4].z])
	else:
		atom1 = array([float(lol[at1][1]), float(lol[at1][2]), float(lol[at1][3])])
		atom2 = array([float(lol[at2][1]), float(lol[at2][2]), float(lol[at2][3])])
		atom3 = array([float(lol[at3][1]), float(lol[at3][2]), float(lol[at3][3])])
		atom4 = array([float(lol[at4][1]), float(lol[at4][2]), float(lol[at4][3])])
	# making appropriate vectors and normals 
	vector1 = atom2-atom1
	vector2 = atom3-atom2
	plane1 = cross(vector1,vector2)
	vector3 = atom2-atom3
	vector4 = atom4-atom3
	plane2 = cross(vector3,vector4)
	# finding dihedral angle
	dihedral = arccos(-dot(plane1,plane2)/(linalg.norm(plane1)*linalg.norm(plane2)))
	# checking the sign of the dihedral then displaying result
	if dot(plane1,vector4) > 0:
		#print degrees(dihedral)
		return  degrees(dihedral)
	else:
		#print -degrees(dihedral)
		return - degrees(dihedral)

### --- get the distance, angle and dihedrals for a Z-matrix --- ###
def getzmat(i):
	line = []
	line.append(getElementName(ifilelol[i+2].e))
	if i > 0:
		line.append(i)
		dist = getdistance(i+1, i+2, ifilelol)
		line.append(dist)
	if i > 1:
		line.append(i-1)
		angle = getangle(i, i+1, i+2, ifilelol)
		line.append(angle)
	if i > 2:
		line.append(i-2)
		dihedral = getdihedral(i-1, i, i+1, i+2, ifilelol)
		line.append(dihedral)
	line.append(-1)
	line.append(-1)
	line.append(-1)
	line.append(-1)
	line.append(-1)
	line.append(-1)
	line.append("\n")
	return line

### --- Get the XYZ coordinates from distance, angle and dihedral data --- ###
def getXYZfromZMAT(lol,at4):
	### Set up the variables to be used for the function
	dist = float(lol[2])
	angle = float(lol[4])
	dihedral = float(lol[6])
	angleRad = radians(angle) # * math.pi / 180
	dihedralRad = radians(dihedral) # * math.pi / 180
	at1 = int(lol[5])-1
	at2 = int(lol[3])-1
	at3 = int(lol[1])-1
	x = 0
	y = 0
	z = 0
	line = []

	### Start to place the atoms in their locations
	if at4 == 0:
		x = 0.00000
		y = 0.00000
		z = 0.00000
	
	elif at4 == 1:
		x = dist
		y = 0.00000
		z = 0.00000
	
	elif at4 == 2:
		a = xyzLOL[at3][1]
		b = dist
		x = a + (dist * cos(math.pi - angleRad))
		y = 0.00000
		z = -dist * sin(math.pi - angleRad)
		
	elif at4 >= 3:
		####The at4 x,y,z coordinates from spherical coord
		Sx = dist * sin(angleRad) * cos(dihedralRad)
		Sy = -dist * sin(angleRad) * sin(dihedralRad) #For some reason I need to have a negative here to get the correct sign in the output..... weird
		Sz = dist * cos(angleRad)
		at4L = [Sx, Sy, Sz]

		###Finding the angle theta
		#Make the list of lists for the three point (z-axis, origin, and translated atom 2) needed for an angle calculation
		Z32 = [[0, 0, 0, 1],  [0, 0, 0, 0],  [0, xyzLOL[at2][1] - xyzLOL[at3][1], xyzLOL[at2][2] - xyzLOL[at3][2], xyzLOL[at2][3] - xyzLOL[at3][3]]]
		#Get theta using the getangle function
		theta = radians(getangle(0, 1, 2, Z32))
		###Rodrigues' rotation formula
		#Create the vectprs needed to calculate k
		vector3 = array([xyzLOL[at3][1], xyzLOL[at3][2], xyzLOL[at3][3]])
		vector2 = array([xyzLOL[at2][1], xyzLOL[at2][2], xyzLOL[at2][3]])
		vector0 = array([0, 0, 1])
		#Calculate k for the Rodrigues rotation formula
		k = cross((vector2-vector3), vector0)/linalg.norm(cross((vector2-vector3), vector0))
		#Generate an array for translated 1
		T1 = [(xyzLOL[at1][1]-xyzLOL[at3][1]), (xyzLOL[at1][2]-xyzLOL[at3][2]), (xyzLOL[at1][3]-xyzLOL[at3][3])]
		#Calculate the Rodrigues rotation matrix
		RR23T1 = dot(T1, cos(theta)) + dot(cross(k,T1), sin(theta)) + dot(dot(k,(dot(k,T1))), (1-cos(theta)))
		#Make the list of lists for the four points (x-axis, z-axis, origin, and rotated translated 1) needed for a dihedral calculation
		XZ31 = [[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, RR23T1[0], RR23T1[1], RR23T1[2]]]
		#Get phi using the getdihedral function
		phi = radians(getdihedral(0,1,2,3,XZ31))
		###Rotation matrix 
		#Create the array for the rotation matrix including dihedral phi
		RM = array([[cos(phi), sin(phi), 0], [-sin(phi), cos(phi), 0], [0, 0, 1]])
		#Calculate the dot product of the rotation matrix and the coordinates for 4 (from spherical)
		RM4 = dot(RM, at4L)
		#Calculate the rotated coordinates of the rotated coordinates of atom 4
		RRN23RM4 = dot(RM4, cos(-theta)) + dot(cross(k,RM4), sin(-theta)) + dot(dot(k,(dot(k,RM4))), (1-cos(-theta)))
		#Final coordinates that are rotated, rotated and translated
		x = RRN23RM4[0] + xyzLOL[at3][1]
		y = RRN23RM4[1] + xyzLOL[at3][2]
		z = RRN23RM4[2] + xyzLOL[at3][3]
	#Putting everything into a list to send back
	line.append(lol[0])
	line.append(x)
	line.append(y)
	line.append(z)
	return line

### --- Open parent file --- ###
f = open(ifile)
ifileList = f.readlines()
f.close()
### --- Create some variables needed for parsing the input file --- ###
# ifilelength is used to gauge what line we are on and also allows us to determine the number of lines in the file
# ifilelol is a list of lists (lol) that will hold all the data of the input file
# a_list is a temp list that will be appended to the ifilelol
################################
ifileLength = len(ifileList)
#iFileAtomNum = ifileLength - 2
ifilelol = [0] * ifileLength
a_list = []
#### --- Input/parse the input file into a list of lists called ifilelol --- ###
for i in range(ifileLength):
	if i == 0:
		ifilelol[i] = int(ifileList[i].rstrip())
		iFileAtomNum = int(ifileList[i].rstrip())
	elif i == 1:
		ifilelol[i] = str(ifileList[i].rstrip())
	elif i >= 2:
		line = ifileList[i].rstrip().split()
		ifilelol[i] = Atom()
		ifilelol[i].e = getElementNum(line[0])
		ifilelol[i].x = float(line[1])
		ifilelol[i].y = float(line[2])
		ifilelol[i].z = float(line[3])

### --- Get bonds --- ###
covBonds = []
covHBonds = []
covTMBonds = []
nearestNeighbor = []
neighborStart = [0, 1000000]
### --- Generate bond lists --- ###
for i in range(0,int(iFileAtomNum)):
	nearestNeighbor.append(list(neighborStart))
	for j in range(0,int(iFileAtomNum)):
		if i != j:
			distij = getdistance(i+2, j+2, ifilelol)
			if j > i:
				if distij <= 2.25 and ifilelol[i+2].e != 1 and ifilelol[j+2].e != 1 and ((getElementName(ifilelol[i+2].e).lower() not in elementLarge) or (getElementName(ifilelol[j+2].e).lower() not in elementLarge)):
				#	distList = [i+1, j+1, distij]
					ifilelol[i+2].neighbors.append(j+1)
					ifilelol[j+2].neighbors.append(i+1)
				#	covBonds.append(distList)
					#print str(i+2) + "\t" + str(j+2) + "\t" + str(distij)
				elif distij <= 3 and ((getElementName(ifilelol[i+2].e).lower() in elementLarge) and (getElementName(ifilelol[j+2].e).lower() in elementLarge)):
				#	distList = [i+1, j+1, distij]
					ifilelol[i+2].neighbors.append(j+1)
					ifilelol[j+2].neighbors.append(i+1)
				#	covTMBonds.append(distList)
				elif distij <= 1.3 and (ifilelol[i+2].e == 1 or ifilelol[j+2].e == 1):
				#	distList = [i+1, j+1, distij]
					ifilelol[i+2].neighbors.append(j+1)
					ifilelol[j+2].neighbors.append(i+1)
				#	covHBonds.append(distList)
			if distij < nearestNeighbor[i][1]:
			#	nearestNeighbor[i][0] = j + 1
			#	nearestNeighbor[i][1] = distij
				ifilelol[i+2].nearest = j+1
### --- Remove hydrogen bonds from bond list --- ###
#for i in range(0,len(covHBonds)):
#	if (covHBonds[i][0] != nearestNeighbor[covHBonds[i][0]][0]) and (covHBonds[i][1] != nearestNeighbor[covHBonds[i][0]][0]):
#		del covHBonds[i]	

#print "Covalent bonds to Hydrogen:"
#print covHBonds
#print "Covalent bonds between \"heavy\" atoms."
#print covBonds
#print "Covalent bonds between TM atoms."
#print covTMBonds
#print "Each Atoms nearest neighbor."
#print nearestNeighbor

### --- Reorder neighbor list to be in numerical order --- ###
for i in range(2,ifileLength):
	for j in range(len(ifilelol[i].neighbors)):
		for k in range(len(ifilelol[i].neighbors)):
			if ifilelol[ifilelol[i].neighbors[k]+1].e < ifilelol[ifilelol[i].neighbors[j]+1].e:
				ifilelol[i].neighbors.insert(0, ifilelol[i].neighbors.pop(k))

### --- Determine if the atom is in a ring --- ### NOT DONE!!!!
def isinring(at1, lol):
	neighborList = []
	for i in range(len(lol[at1].neighbors)):
		#print lol[lol[at1].neighbors[i]+1].neighbors
		neighborList.append(lol[lol[at1].neighbors[i]+1].neighbors)
	#print neighborList
	### --- Search for rings --- ###	
	ringSearchDone = "no"
	while ringSearchDone == "no":
		for j in range(len(neighborList)):
			for i in range(len(neighborList[j])):
				print lol[neighborList[j][i]+1].neighbors
				#neighborList.append(lol[lol[neighborList[j][i]+1].neighbors[i]+1].neighbors)
		print neighborList
		ringSearchDone = "yes"
	return

### --- Output strings to discern certain types of atoms based on their neighbors --- ###
ffsmilelol = []
for i in range(2,ifileLength):
	line = [str(getElementName(ifilelol[i].e))]
	###line = [getElementName(ifilelol[i].e), ifilelol[i].x, ifilelol[i].y, ifilelol[i].z, ifilelol[i].neighbors, ifilelol[i].neighborsdist, ifilelol[i].nearest, ifilelol[i].charge]
	for atom in ifilelol[i].neighbors:
		line.append(" | ")
		line.append(str(getElementName(ifilelol[atom+1].e)))
		for subatoms in ifilelol[atom+1].neighbors:	
			if subatoms + 1 != atom:
				line.append("(")
				line.append(str(getElementName(ifilelol[subatoms+1].e)))
				if len(ifilelol[subatoms+1].neighbors) > 1:
					line.append("(")	
					for subsubatoms in ifilelol[subatoms+1].neighbors:
						if subsubatoms + 1 != subatoms:
							line.append(str(getElementName(ifilelol[subsubatoms+1].e)))
					line.append(")")
				line.append(")")
	ffsmilelol.append(''.join(str(e) for e in line).split('|'))
#	print ''.join(str(e) for e in line)
	line[:] = []

print ffsmilelol
#isinring(2, ifilelol)



#### --- Output file in Z-matrix format --- ###
#f = open(ofile, 'w+')
#f.write(str(iFileAtomNum) + "\n")
#f.write(ifilelol[1] + "\n")
#for i in range(0,int(iFileAtomNum)):
#	linetemp = [x for x in zmatlol[i] if x != -1]
#	line = '\t'.join(str(e) for e in linetemp) #, e != -1)
#	f.write(line)
#f.close()


#C   0.0   0.0   0.0
#H   1.08667477625   0.0   0.0
#C   -0.710304820313   0.0   -1.19862882797
#H   -0.190303313054   0.00285187546684   -2.15257116105
#C   -2.10196298123   0.00218056792464   -1.16393829937
#H   -2.63732817294   0.0171920874033   -2.10597724946
#C   -2.81964532451   -0.00437264458075   0.0480714706799
#P   -4.6702532026   -0.0170742413099   0.00664801286375
#Pd   -5.54865169312   0.128925626685   2.17205850687

######################################################################
### END OF SCRIPT	
######################################################################