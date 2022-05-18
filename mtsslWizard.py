"""
*************************************************************
*	This is a command line version of mtsslWizard			*
*	Example: python mtsslWizard.py 2lzm.pdb A [89,85] r1	*
*	Version 0.1 alpha, 30.10.2018							*
*	Gregor Hagelueken										*
*	hagelueken@uni-bonn.de			     					*
*************************************************************
"""
from mWclasses.label import Label as MWlabel
from mWclasses.ensemble import Ensemble as MWensemble
from mWclasses.distanceDistribution import DistanceDistribution as MWdistanceDistribution
from mWclasses.distanceMap import DistanceMap as MWdistanceMap
from Bio.PDB import *
import numpy
import itertools
import os
import sys
import ast
import random
import scipy.spatial.distance
import copy

def labelMany(pdbFile, distances):
	for distance in distances:
		for labelPosition in distance:
			residueId = labelPosition["resi"]
			chainId = labelPosition["chain"]
			label = labelPosition["label"]
			numberOfTries = 100000
			totalNumberOfTries = 100000
			triesSoFar = 0
			parser = PDBParser()
			proteinA = parser.get_structure("proteinA", pdbFile)
			proteinAcoords = []
			for model in proteinA:
				for chain in model:
					for residue in chain:
						if not residue.id == 'W' or not 'H_' in residue.id: 
							if not (residue.id[1] == residueId and chain.id == chainId):
								for atom in residue:
									proteinAcoords.append(atom.get_coord())
			proteinAcoords = numpy.asarray(proteinAcoords)
			mWlabel = MWlabel.fromfile("labels/%s.txt" %label)
			labelStructure = parser.get_structure("label", os.path.dirname(os.path.abspath(__file__)) + "/mWclasses/labels/%s.pdb" %label)
	
			labelAtoms = labelStructure.get_atoms()
			labelCoords = []
			labelSuper = []
			for atomName in mWlabel.atomNames:
				labelAtoms = labelStructure.get_atoms()
				for atom in labelAtoms:
					if atom.get_name() == atomName:
						labelCoords.append(atom.get_coord())
	
			for atomName in mWlabel.atomsForSuperposition:
				labelAtoms = labelStructure.get_atoms()
				for atom in labelAtoms:
					if atom.get_name() == atomName:
						labelSuper.append(atom)

			proteinSuper = []
			chain = proteinA[0][chainId]

			residue = chain[int(residueId)]
			for labelAtom in labelSuper:
				for residueAtom in residue:
					if residueAtom.get_name() == labelAtom.get_name():
						proteinSuper.append(residueAtom)
	
			super_imposer = Superimposer()
			super_imposer.set_atoms(proteinSuper, labelSuper)
			super_imposer.apply(labelStructure.get_atoms())
			print "	   r.m.s. of label superposition: %1.3f Ang." %super_imposer.rms
	
			#get new coords after superposition, put in correct order
			labelAtoms = labelStructure.get_atoms()
			labelCoords = []

			for atomName in mWlabel.atomNames:
				labelAtoms = labelStructure.get_atoms()
				for atom in labelAtoms:
					if atom.get_name() == atomName:
						labelCoords.append(atom.get_coord())

			mWlabel.movingAtoms = labelCoords
			mWensemble = MWensemble()
			mWensemble.name = "mW"
			numberOfProcesses = 4
			rotamers = mWlabel.generateEnsembleMulti(mWlabel.movingAtoms, proteinAcoords, mWlabel.numberToFind['thorough'], numberOfTries, totalNumberOfTries, triesSoFar, mWensemble.name, 2.5, 5)
			triesSoFar += numberOfTries
			mWensemble.rotamers = rotamers
			spinPositions = []
			for rotamer in rotamers:
				for idx, atomName in enumerate(mWlabel.atomNames):
					if atomName == mWlabel.spinLocation:
						coordinate = rotamer.atoms[atomName].coordinate
						spinPositions.append(coordinate)
			labelPosition["spinPositions"] = spinPositions
			numberOfRotamers = len(rotamers)
			print "	   found %i rotamers" %len(rotamers)
			filename = "%s_%s_%s.pdb" %(pdbFile, chainId, residueId)
			mWensemble.writePDB(filename = filename)
	print "Done!"
	return distances

def label(pdbFile, chainId, residueIds, labelname):
	for residueId in residueIds:
		print "Labelling %s chain: %s residue: %i with: %s" %(pdbFile, chainId, residueId, labelname)
		numberOfTries = 100000
		totalNumberOfTries = 100000
		triesSoFar = 0
		parser = PDBParser()
		proteinA = parser.get_structure("proteinA", pdbFile)
		proteinAcoords = []
		for model in proteinA:
			for chain in model:
				for residue in chain:
					if not residue.id == 'W' or not 'H_' in residue.id: 
						if not (residue.id[1] == residueId and chain.id == chainId):
							for atom in residue:
								proteinAcoords.append(atom.get_coord())
		proteinAcoords = numpy.asarray(proteinAcoords)
		mWlabel = MWlabel.fromfile("labels/%s.txt" %labelname)
		labelStructure = parser.get_structure("label", "mWclasses/labels/%s.pdb" %labelname)
	
		labelAtoms = labelStructure.get_atoms()
		labelCoords = []
		labelSuper = []
		for atomName in mWlabel.atomNames:
			labelAtoms = labelStructure.get_atoms()
			for atom in labelAtoms:
				if atom.get_name() == atomName:
					labelCoords.append(atom.get_coord())
	
		for atomName in mWlabel.atomsForSuperposition:
			labelAtoms = labelStructure.get_atoms()
			for atom in labelAtoms:
				if atom.get_name() == atomName:
					labelSuper.append(atom)

		proteinSuper = []
		chain = proteinA[0][chainId]

		residue = chain[int(residueId)]
		for labelAtom in labelSuper:
			for residueAtom in residue:
				if residueAtom.get_name() == labelAtom.get_name():
					proteinSuper.append(residueAtom)
	
		super_imposer = Superimposer()
		super_imposer.set_atoms(proteinSuper, labelSuper)
		super_imposer.apply(labelStructure.get_atoms())
		print "	   r.m.s. of label superposition: %1.3f Ang." %super_imposer.rms
	
		#get new coords after superposition, put in correct order
		labelAtoms = labelStructure.get_atoms()
		labelCoords = []

		for atomName in mWlabel.atomNames:
			labelAtoms = labelStructure.get_atoms()
			for atom in labelAtoms:
				if atom.get_name() == atomName:
					labelCoords.append(atom.get_coord())

		mWlabel.movingAtoms = labelCoords
		mWensemble = MWensemble()
		mWensemble.name = "mW"
		numberOfProcesses = 1
		rotamers = mWlabel.generateEnsembleMulti(mWlabel.movingAtoms, proteinAcoords, mWlabel.numberToFind['thorough'], numberOfTries, totalNumberOfTries, triesSoFar, mWensemble.name, 3.4, 0)
		triesSoFar += numberOfTries
		mWensemble.rotamers = rotamers
		numberOfRotamers = len(rotamers)
		print "	   found %i rotamers" %len(rotamers)
		filename = "%s_%s_%s.pdb" %(pdbFile, chainId, residueId)
		mWensemble.writePDB(filename = filename)
	print "Done!"

def quick_map(atoms1, atoms2):
	dist = scipy.spatial.distance.cdist(atoms1, atoms2)
	return dist.flatten()

print sys.argv
if len(sys.argv) == 5:
	pdbFile = sys.argv[1]
	chainId = sys.argv[2]
	residueIds = ast.literal_eval(sys.argv[3])
	labelname = sys.argv[4]
	print "parameters: ", pdbFile, chainId, residueIds, labelname
	try:
		label(pdbFile, chainId, residueIds, labelname)
	except:
		print "something went wrong!"
elif len(sys.argv) == 2:
	normalMode = False
	if normalMode:
		inputFile = sys.argv[1]
		with open(inputFile) as f:
			content = f.readlines()
		pdbFiles = []
	
		for line in content:
			if ".pdb" in line:
				pdbFiles.append(line.rstrip())
		content = content[len(pdbFiles):]
		for pdbFile in pdbFiles:
			distances = []
			for line in content:
				labelPositions = []
		
				thisLine = line.split(",")
				for item in thisLine:
					labelPosition = {}
					resiAndChain = item.split("_")
					labelPosition["chain"] = resiAndChain[0].rstrip()
					labelPosition["resi"] = resiAndChain[1].rstrip()
					labelPosition["label"] = resiAndChain[2].rstrip()
					labelPositions.append(labelPosition)
				distances.append(labelPositions)
			distances = labelMany(pdbFile, distances)
			for distance in distances:
				#get ensembles
				dist = []
				ensembleAtoms = []
				for ensemble in distance:
					coordinates = []
					coordinates = ensemble["spinPositions"]
					ensembleAtoms.append(numpy.asarray(coordinates))
				for pair in itertools.combinations(ensembleAtoms,2):
					mwDistanceDistribution = MWdistanceDistribution()
					pair_distances = mwDistanceDistribution.calculateDistanceDistribution(pair[0], pair[1])
					dist.extend(pair_distances)
				histogram = numpy.histogram(dist, numpy.arange(100))
				envelopePlot = numpy.zeros((100,2))
				envelopePlot[0:99] = numpy.column_stack((histogram[1][0:len(histogram[1])-1], histogram[0]))
	
				#put point in mid of bin
				envelopePlot[:,0] += 0.5 
				normEnvelopePlot = numpy.copy(envelopePlot)
				normEnvelopePlot[:,1] = normEnvelopePlot[:,1]/numpy.amax(histogram[0])
	
				#combine dist and histogram to single array before output
				output = numpy.column_stack((envelopePlot, normEnvelopePlot[:,1]))
				averageDistance = numpy.average(dist)
				print averageDistance
				distributionString = "["
				for row in output:
					#print row
					x = row[0]
					y = row[2]
					#print x, y
					newPoint = "{x:%1.2f, y:%1.2f}," %(x, y)
					distributionString += newPoint
				distributionString += "]"
				csvString = "%s\n" %averageDistance
				for row in output:
					#print row
					x = row[0]
					y = row[2]
					#print x, y
					newPoint = "%1.2f\t%1.2f\n" %(x, y)
					csvString += newPoint
				labelSites = ""
				for ensemble in distance:
					labelSites += ensemble["chain"]
					labelSites += ensemble["resi"]
					labelSites += ensemble["label"]
					labelSites += "_"
				filename = "%s-%s.txt" %(pdbFile, labelSites)
				with open(filename, "w") as text_file:
					text_file.write(csvString)
	else:
		inputFile = sys.argv[1]
		with open(inputFile) as f:
			content = f.readlines()
		pdbFiles = []
	
		for line in content:
			if ".pdb" in line:
				pdbFiles.append(line.rstrip())
		content = content[len(pdbFiles):]
		for pdbFile in pdbFiles:
			distances = []
			for line in content:
				labelPositions = []
		
				thisLine = line.split(",")
				for item in thisLine:
					labelPosition = {}
					resiAndChain = item.split("_")
					labelPosition["chain"] = resiAndChain[0].rstrip()
					labelPosition["resi"] = resiAndChain[1].rstrip()
					labelPosition["label"] = resiAndChain[2].rstrip()
					labelPositions.append(labelPosition)
				distances.append(labelPositions)
			requestedDistances = labelMany(pdbFile, distances)
			weightsOfFirstAtom = [0,10,50,100]
			#numberOfTrialAtoms = 1000
			for weightOfFirstAtom in weightsOfFirstAtom:
				avgErrors = []
				for j in range (1000): #100 "experiments"
					#requestedDistances is a list of the requested distributions from the input file.
					#choose a random distance
					thisDistance = random.choice(requestedDistances)
				
					#get the spin Positions
					trialAtoms1 = copy.copy(thisDistance[0]["spinPositions"])
					trialAtoms2 = copy.copy(thisDistance[1]["spinPositions"])
					
					avgCoordTrialAtoms1 = numpy.mean(trialAtoms1, axis = 0)
					avgCoordTrialAtoms2 = numpy.mean(trialAtoms2, axis = 0)
					avgDist = numpy.linalg.norm(avgCoordTrialAtoms1 - avgCoordTrialAtoms2)
					fixedPosition1 = random.choice(trialAtoms1)
					fixedPosition2 = random.choice(trialAtoms2)
					for i in range(1, int(weightOfFirstAtom/100.0 * len(trialAtoms1))):
						trialAtoms1[i] = fixedPosition1
					
					for i in range(1, int(weightOfFirstAtom/100.0 * len(trialAtoms2))):
						trialAtoms2[i] = fixedPosition2

					distances = quick_map(trialAtoms1, trialAtoms2)
					errors = distances - avgDist
					avgErrors.append(numpy.mean(errors))
				numpy.savetxt("errors_%i.txt" %(weightOfFirstAtom), avgErrors)

else:
	print "Example: python mtsslWizard.py 2lzm.pdb A [89,85] r1"
