#!/usr/bin/python
import sys
import commands

################
#
#This is a code to interpret the output from a CRYSTAL09 calculation into a usable format to interface with the NBO code of BDD
#The output of this will still have to be processed by a FORTRAN code in order to make it compatible with NBO
#That code is called process_crystal.exe (source=crystal_nbo.f90)
#
####BDD 6/22/12#


fn = sys.argv[1]


##############################
#This function is meant to take a single integer from a file
#First it grabs the line on which 'keyword' is located. 
#Then shift is how many 'words' from the end the paramter is on the line. Shift=1 corresponds to the last element.
def get_param(file_name,keyword,shift):
	call = "grep '"  + keyword + "' " + file_name  #prepare the grep call for 'keyword'
	line = commands.getoutput( call )
	words = line.split()   #break the line into words (by spaces) so that only the value can be obtained
	return int(words[len(words)-shift])   #Convert the word to an integer
##############################


############Start by finding out some general information about the system###############
natom = get_param(fn, 'ATOMS IN THE UNIT CELL', 1) #Number of atoms in the unit cell
nkpts = get_param(fn, 'K POINTS IN THE IBZ', 1) #Number of k-points used 
nbasis = get_param(fn, 'NUMBER OF AO', 7) #Numer of basis functions

#Output all dem parameters
print nkpts, '  #nkpts'
print natom, '  #natom'
print nbasis, '  #nbasis'
print


###############################
#This functions will look for a keyword 'test' in a file and then collect 'dim' number of lines into 'arr' after skipping 'shift' lines
def parse(file_name,test,shift,dim,arr):
	file_name.seek(0)  #Start by shifting the cursor in the file to the beginning, this way guaranteed to find the keyword
	j=0    #Counter that will be used to know when to start looking at lines
	for line in file_name:
		if test in line:
			j=dim+shift #When the keyword has been found, set up the counter to know when to start and how many lines to print
		elif j > 0:  #Start iterating counter
			if j <= dim:  #if the shift interpreting lines from the file
				arr[dim-j]=line.split()  #Lines are read in as split into words, greatly easing manipluation outside of the parse functions
			j = j-1
###############################


######Now we will begin parsing out information besides systems parameters#######

#Begin by casting the CRYSTAL output file as a file type
fn = file( sys.argv[1], 'r')


##First thing to find is the lattice vectors##
lattice=range(3)
parse(fn,'CARTESIAN COMPONENTS',1,3,lattice)
#print 'lattice vectors in ang, each line is new vector, colums are cartesian dimension'
#for i in range(3):  #Write out the lattice vectors to the screen
#	for j in range(3):
#		print lattice[i][j],
#	print
#print


##Now find some information on the atoms of the system##
temp=range(natom) #The line has various useful information so this array will be further interpreted
parse(fn,'CARTESIAN COORDINATES',3,natom,temp)

#Set up the arrays to store all the different information in coord_temp
atom_number = range(natom)  #array for atomic numbers
atom_sym = range(natom)     #array for atomic symbold
atom_coord = range(natom)  #array for atomic numbers

#Actually get the information out of coord_temp
for x in range(natom):
 	#print x, coord_temp[x][1]
	atom_number[x] = temp[x][1]
	atom_sym[x] = temp[x][2]
	atom_coord[x] = temp[x][3:6]


##Get information about the basis set##
fn.seek(0)

for line in fn:
	if 'END' in line:  #made it to beginning of basis set block
		break
temp=range(0) #This will hold all basis set information
for line in fn:
	if 'END' in line:
		break
	else:
		temp.append(line.split())

#print 'processed basis input'
#print temp
#print

x = 0
itype = 0

shells = range(itype)
dim = range(itype)
index = range(itype)
orb_types = range(itype)
num_gauss = range(itype)

shift = 0

while x < 1:

	index.append(int(temp[shift][0]))
	shells.append(int(temp[shift][1]))
	dim.append(0)

        orb_types.append(range(shells[itype]))
	num_gauss.append(range(shells[itype]))

	#print 'itype '+str(itype)+':',index[itype],shells[itype] #,dim[itype]

	for j in range(shells[itype]):
		shift = shift + 1

		orb_types[itype][j] = int(temp[shift][1]) #Type of basis function

		if orb_types[itype][j] == 0:
			dim[itype] = dim[itype] + 1  #s type
                elif orb_types[itype][j] == 1:
                        dim[itype] = dim[itype] + 4  #sp type
                elif orb_types[itype][j] == 2:
                        dim[itype] = dim[itype] + 3  #p type
                elif orb_types[itype][j] == 3:
                        dim[itype] = dim[itype] + 5  #d type
		else:
			sys.exit('Unrecognized orbital index in basis set input '+orb_types[itype][1])

		num_gauss[itype][j] = int(temp[shift][2])
		for k in range(num_gauss[itype][j]):
			shift = shift + 1

	shift = shift + 1 #move to next potential type

	if temp[shift][0] == '99' and temp[shift][1] == '0': #This is must be at the end of the CRYSTAL basis set entry, is always last line
		x = 1  #the while loop is now exited
	else:
		itype = itype + 1

#print

ntype = len(index)	

#print 'index     ',index
#print 'shells    ',shells
#print 'orb_type  ',orb_types
#print 'num_gauss ',num_gauss
#print 'ntype     ',ntype
#print 'dim       ',dim
#print


#Should have everything I need out of the file by here so close it
fn.close()
#print fn.closed

###########################
def clean_str(ls):
        """Create a string of all the elements in a list, without commas or brackets"""
        cln = ''
        for x in range(len(ls)):
                if x < len(ls)-1:
                        cln = cln + str(ls[x]) + '  '
                else:
                        cln = cln + str(ls[x])   #Avoid a trailing space at the end of the string
        return cln
###########################



#Now take the exponen and coeff out of temp, into a controlled array
basis_output = range(ntype)
shift = 0
for i in range(ntype):
	basis_output[i] = range(shells[i])
	for j in range(shells[i]):
		shift = shift + 1
		basis_output[i][j] = range(num_gauss[i][j])
		for k in range(num_gauss[i][j]):
			shift = shift + 1
			basis_output[i][j][k] = temp[shift]
		#print
	shift = shift + 1



#Different basis functions can be used for atoms of the same element type in CRYSTAL
#This code assumes that is not the case and this loop will exit the code if that is so
if ntype > 1:
	for i in range(ntype):
		for j in range(i+1,ntype):
			if index[i]%100 == index[j]%100:
				sys.exit('code is not set up to support different basis sets for same element type')


#Make sure the indices calculated here cover all of the atomic numbers gathered above
atom_types=range(natom)
dim_test = 0
for i in range(natom):
	for j in range(ntype):
		if index[j] == int(atom_number[i]):
			atom_types[i]=j
			dim_test = dim_test + dim[j]
			break
	else:
		sys.exit('For atom '+str(i)+' there is no corresponding basis set '+atom_number[i])


#Make sure the same number of basis functions was calcualted out as was taken from output file
if dim_test != nbasis:
	sys.exit('The total numer of basis functions '+str(dim_test)+' is not the same as that taken from file '+str(nbasis))


#Make sure all atom types in basis set are represented in the system.
for j in range(ntype):
	for i in range(natom):
		if index[j] == int(atom_number[i]):
			break
	else:
		sys.exit('Type '+str(j)+' in the basis set does not correspond to any atom.')



ibasismap=range(0)
shift=1

ishell = 0
ishellmap=range(0)

for i in range(natom):
	itype = atom_types[i]
	ibasismap.append(shift)
	for j in range(shells[itype]):

		ishell = ishell + 1

		sp = 0

		if orb_types[itype][j] == 0:  #s-type
			l = 0
			shell_num = 1
		elif orb_types[itype][j] == 1:  #sp-type

                        print num_gauss[itype][j],'   num gauss'
                        for x in range(num_gauss[itype][j]):
                                print basis_output[itype][j][x][0],' ',
                        print '   alphas'
                        for x in range(num_gauss[itype][j]):
                                print basis_output[itype][j][x][1],' ',
                        print '   coeffs'
                        print clean_str(atom_coord[i]), ' position'
			print i+1, '   atom'
                        print 0, 1, '    #l and m'
                        print

			l = 1
			sp = 1
			shell_num = 4
		elif orb_types[itype][j] == 2:  #p-type
			l = 1
			shell_num = 3
		elif orb_types[itype][j] == 3:  #d-type
			l = 2
			shell_num = 5
		for m in range(2*l + 1):

			print num_gauss[itype][j],'   num gauss'
			for x in range(num_gauss[itype][j]):
				print basis_output[itype][j][x][0],' ',
			print '   alphas'
			for x in range(num_gauss[itype][j]):
				print basis_output[itype][j][x][1+sp],' ',
			print '   coeffs'
			print clean_str(atom_coord[i]),' position'
			print i+1, ' atom'
			print l, m+1, '    #l and m'
			print


		shift = shift + shell_num
		for x in range(shell_num):
			ishellmap.append(ishell)

print

ibasismap.append(shift)



#Write out the atomic information to the screen
print clean_str(atom_number),'    #iatnum'
print clean_str(atom_sym),'    #symbols'
print clean_str(ibasismap),'    #ibasismap'
print clean_str(ishellmap),'    #ishellmap'
print

#print 'lattice vectors in ang, each line is new vector, colums are cartesian dimension'
for i in range(3):  #Write out the lattice vectors to the screen
	print clean_str(lattice[i])
        #for j in range(3):
        #        print lattice[i][j],
        #print
print


for x in range(natom):
	print clean_str(atom_coord[x])
        #for y in range(3):
        #        print atom_coord[x][y],
        #print
print


