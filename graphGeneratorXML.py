################################### Graph Generator v1.0 ################################
#											#
# 				    University of York					#
# 				     Graceful Project					#
# 				  Created by Pedro Campos				#
# 					York, 2015					#
#											#
#########################################################################################

import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import re
import os, errno, sys
import random
import math
import string
import cPickle as pickle

#	OBJECTS USED IN THE GRAPH PARSING AND TRAFFIC GENERATION	#

class Node:
	def __init__(self, nodeId, name):
		self.nodeId = nodeId
		self.name = name

class Edge:
	def __init__(self, sourceId, targetId, networkLoad):
		self.sourceId = sourceId
		self.targetId = targetId
		self.networkLoad = networkLoad

class TrafficEdge:
	def __init__(self, sourceId, targetId, triggerPacketSize):
		self.sourceId = sourceId
		self.targetId = targetId
		self.triggerPacketSize = triggerPacketSize

class TrafficTarget:
	def __init__(self, nodeId, outputPackets, targetId, dataSourcesList, triggerPacketsList):
		self.nodeId = nodeId
		self.outputPackets = outputPackets
		self.targetId = targetId
		self.dataSourcesList = dataSourcesList
		self.triggerPacketsList = triggerPacketsList

class TrafficElement:
	def __init__(self, nodeId, targetsList):
		self.nodeId = nodeId
		self.targetsList = targetsList

class TrafficList:
	def __init__(self, graphID, trafficElements, maxDegree, networkBandwidths):
		self.graphID = graphID
		self.trafficElements = trafficElements
		self.maxDegree = maxDegree
		self.networkBandwidths = networkBandwidths

class basicGraph:
	def __init__(self, ID, nodeList, edgeList):
		self.ID = ID
		self.nodeList = nodeList
		self.edgeList = edgeList

class Graph:
	def __init__(self, ID, nodeList, edgeList, trafficList):
		self.ID = ID
		self.nodeList = nodeList
		self.edgeList = edgeList
		self.trafficList = trafficList



# Function which imports graphs written in XML and parses them to objects, returning a basicGraph object
def importXML(filename, graphID=1):

	nodeList = []
	edgeList = []

	activeNodeList = False
	activeEdgeList = False
	activeNode = False
	activeEdge = False

	fOpen = open(filename, 'r')
	data = fOpen.readlines()
	# Dump everything into list, so XML can be added later.
	dumpData=data
	nodeID = 0

	for line in data:
		
		# Lines of interest #

		nListOpen = re.findall(r'<mc:NodeList>.*', line)
		nListClose = re.findall(r'</mc:NodeList>.*', line)
		eListOpen = re.findall(r'<mc:EdgeList>.*', line)
		eListClose = re.findall(r'</mc:EdgeList>.*', line)
		nodeOpen = re.findall(r'<mc:Node>.*', line)
		nodeClose = re.findall(r'</mc:Node>.*', line)
		edgeOpen = re.findall(r'<mc:Edge>.*', line)
		edgeClose = re.findall(r'</mc:Edge>.*', line)

		nodeIdLine = re.findall(r'<mc:id>.*', line)
		nodeNameLine = re.findall(r'<mc:name>.*', line)
		edgeSourceLine = re.findall(r'<mc:sourceId>*', line)
		edgeTargetLine = re.findall(r'<mc:targetId>.*', line)
		edgeNetworkLine = re.findall(r'networkLoad.*', line)

		numbers = re.findall(r'\d+', line)

		# Activation lines #

		if nListOpen:
			activeNodeList = True
		if nListClose:
			activeNodeList = False
		if eListOpen:
			activeEdgeList = True
		if eListClose:
			activeEdgeList = False
		if nodeOpen:
			activeNode = True
		if nodeClose:
			activeNode = False
		if edgeOpen:
			activeEdge = True
		if edgeClose:
			activeEdge = False
		
		# XML parsing #

		if activeNodeList and activeNode and nodeIdLine:
			nodeID = repr(int(numbers[0])-1)
		if activeNodeList and activeNode and nodeNameLine:
			l = line.split('<mc:name>')
			nodeName = l[1].split('</mc:name>')
			nodeList.append(Node(nodeID, nodeName[0]))

		if activeEdgeList and activeEdge and edgeSourceLine:
			sourceId = repr(int(numbers[0])-1)
		if activeEdgeList and activeEdge and edgeTargetLine:
			targetId = repr(int(numbers[0])-1)
		if activeEdgeList and activeEdge and edgeNetworkLine:
			networkLoad = numbers[0]
			edgeList.append(Edge(sourceId, targetId, networkLoad))


	return basicGraph(graphID, nodeList, edgeList), dumpData
		

# function which takes a graph in the form of objects and returns a TrafficList object
def addRandomTraffic(graph, graphID, connectivity=0.5, dependenceConstant=0, localThreshold=0.5, outputSizelocalThreshold=0.5, inputPacketRange=1, outputPacketRange=1, maxDegree=5):
	
	sampleSize = 2000
	
	# Generate the normal distribution from which numbers get picked
	outDist = np.random.normal(loc=1, scale=math.sqrt(outputPacketRange), size=sampleSize)
	inDist = np.random.normal(loc=1, scale=math.sqrt(inputPacketRange), size=sampleSize)

	# Round and typecast distribution to get only values from [1, packetRanges]

	castOut = []
	castIn = []

	for i in range(len(outDist)):
		val = int(round(outDist[i]))
		if val > 0:
			castOut.append(val)
		val = int(round(outDist[i]))
		if val > 0:
			castIn.append(val)

	trafficTargets = []	

	for i in graph.nodeList:
		targets = []
		inputs = []
		outputSizeList = []
		for k in graph.edgeList:
			if k.sourceId == i.nodeId:
				targets.append(k.targetId)
			if k.targetId == i.nodeId:
				inputs.append(k.sourceId)
		if not inputs:
			inputs.append('SOURCE')
			#inputs.append('LOCAL')
		if not targets:
			targets.append('SINK')

		connectCheck = [False]*len(inputs)

		trafficEdges = []
		
		# Go through every output and do a first sweep of connecting inputs, later check for the ones left behind		
		for k in range(len(targets)):
			outputConnected = False
			# Defining output packet size for this edge
			outRandom = castOut[int(random.random() * len(castOut))] + dependenceConstant
			outputSizeList.append(outRandom)

			# Defining which inputs are connected to this output		
			while(outputConnected==False):
				for m in range(len(inputs)):
					if random.random() > connectivity:
						#triggerPacketSize = int(random.choice(castIn)) + dependenceConstant
						triggerPacketSize = castIn[int(random.random() * len(castIn))] + dependenceConstant
						connectCheck[m] = True
						trafficEdges.append(TrafficEdge(inputs[m], targets[k], triggerPacketSize))
						outputConnected = True


		# Check for unconnected inputs
		for m in range(len(connectCheck)):
			if (connectCheck[m] == False and (inputs[m] != 'LOCAL')):
				outRandom = int(random.random() * len(targets))
				triggerPacketSize = castIn[int(random.random() * len(castIn))] + dependenceConstant
				connectCheck[m] = True
				trafficEdges.append(TrafficEdge(inputs[m], targets[outRandom], triggerPacketSize))

		
		# LOCAL channel is not currently being used, and is therefore commented out
		# Add randomised local channel
		#if (random.random() > localThreshold and not (any(j=='SOURCE' for j in inputs))):
		#		outRandom = int(random.random() * len(targets))
		#		triggerPacketSize = int(random.choice(castIn)) + dependenceConstant
		#		trafficEdges.append(TrafficEdge('LOCAL', targets[outRandom], triggerPacketSize))


		if any(j == False for j in connectCheck):
			print "Something went wrong.\n"
			sys.exit(0)

		for j in range(len(targets)):
			dataSourcesList = []
			triggerPacketsList = []
			for k in trafficEdges:
				if k.targetId == targets[j]:
					dataSourcesList.append(k.sourceId)
					triggerPacketsList.append(k.triggerPacketSize)
			
			trafficTargets.append(TrafficTarget(i.nodeId, outputSizeList[j], targets[j], dataSourcesList, triggerPacketsList))


	trafficTargets = sorted(trafficTargets, key=lambda TrafficTarget: int(TrafficTarget.nodeId))

	groupedTraffic = []
	groupedTargets = []

	currId = trafficTargets[0].nodeId
	for i in trafficTargets:
		if i.nodeId == currId:
			groupedTargets.append(i)
		else:
			x = TrafficElement(currId, groupedTargets)
			groupedTraffic.append(x)
			groupedTargets = []
			currId = i.nodeId
			groupedTargets.append(i)

	# Deal with last one
	x = TrafficElement(nodeId=currId, targetsList=groupedTargets)
	groupedTraffic.append(x)

	# Print stats as we go along	

	#for x in groupedTraffic:
		#print 'Node ID: ', x.nodeId
		#for i in x.targetsList:
		#	print 'Target ID: ', i.targetId, 'Output Packet Size: ', i.outputPackets, 'Sources: ', i.dataSourcesList, 'Trigger Packet Size: ', i.triggerPacketsList

	tList = TrafficList(graphID=graphID, trafficElements=groupedTraffic, maxDegree=maxDegree, networkBandwidths = [])

	pickleName = "Traffic/dag{0}.data".format(repr(graphID).zfill(4))
	f = open(pickleName, 'w')
	pickle.dump(tList, f)
	f.close()

	return tList


# Function which takes a TrafficList object and writes to an XML file.
def exportXML(filename, graph, dataDump, gvName):

	#seed = random.seed(100)

	# Lines for .gv file, to draw directed graph

	gvLines = []

	#if not os.path.exists('TrafficData'):
#		os.makedirs('TrafficData')

	#gfname		= 'TrafficData/tf_'+filename[:-4] + '.xml'
	#filename 	= 'traffic_'+filename[:-4] + '.xml'
	
	# Expression dictionary
	listBegin 	= '  <TrafficList>\n\n'
	listEnd   	= '  </TrafficList>\n\n'
	elementBegin 	= '    <TrafficElement>\n'
	elementEnd 	= '    </TrafficElement>\n\n'
	targetBegin 	= '      <target>\n'
	targetEnd 	= '      </target>\n\n'
	dependencyBegin = '        <dependency>\n'
	dependencyEnd   = '        </dependency>\n'

	cname = filename[:-4] + '.h'

	f = open(filename, 'w')
	c = open(cname, 'w')

	#for i in range(len(dataDump)):
	#	endGraph = re.findall(r'</Graph>', dataDump[i])

	#	if endGraph:
	#		#dataDump.pop(i)
	#		x = 1
	#	else:
	#		g.write(dataDump[i])

	maxDegree = 0
	for element in graph.trafficElements:

		target_map = []
		target_count = []
		src_map = []
		src_count = [[0]]
		target_index = 0

		for target in element.targetsList:


			if (target.targetId=='SINK'):
				target_map.append(target.targetId)
			else:
				target_map.append(target.targetId)

			target_count.append((target.outputPackets))

			for dependency in range(len(target.dataSourcesList)):
		
				if (target.dataSourcesList[dependency]=='LOCAL' or target.dataSourcesList[dependency]=='SOURCE'):
					if not (target.dataSourcesList[dependency] in src_map):		# Adding SOURCE and LOCAL if not they're not already there
						src_map.append(target.dataSourcesList[dependency])
				else:
					if not (target.dataSourcesList[dependency] in src_map):
						src_map.append(target.dataSourcesList[dependency])

				#src_count[target_index].append(target.triggerPacketsList[dependency])

			target_index += 1
			#src_count.append([])


		target_index = 0
		src_map.sort()
		
		for target in element.targetsList:
			srcCts = [-1]*len(src_map)
			for dependency in range(len(target.dataSourcesList)):
				if (target.dataSourcesList[dependency] in src_map):
					srcCts[src_map.index(target.dataSourcesList[dependency])] = target.triggerPacketsList[dependency]
					#if graph.graphID == 20 and element.nodeId == '4':
					#	print target.triggerPacketsList[dependency], src_map.index(target.dataSourcesList[dependency])
					#	print srcCts
			src_count[target_index] = srcCts
			target_index += 1
			src_count.append([])

		src_count.pop(-1)

		if len(target_map) > maxDegree:
			maxDegree = len(target_map)
		if len(src_map) > maxDegree:
			#print src_map
			maxDegree = len(src_map)

	print "Max Degree: ", maxDegree, " for file", filename
	
	c.write("\t//%s\n\n" % (filename[-11:-4]))
	c.write("\tint map[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};\n\n")
	c.write("\t#define MAX_DEGREE %d\n\n" % (maxDegree))
	c.write("\tint t_info[] =\n\t{\n\n")

	numberNodes = 0

	networkCoefficients = []

	# Begin writing traffic XML
	f.write(listBegin)

	sourceNodes = []

	for element in graph.trafficElements:

		numberNodes += 1

		target_map = []
		target_count = []
		src_map = []
		src_count = [[]]
		target_index = 0

		f.write(elementBegin)
		f.write('      <nodeId>{0}</nodeId>\n\n'.format(element.nodeId))
		for target in element.targetsList:

			f.write(targetBegin)
			f.write('        <targetId>{0}</targetId>\n'.format(target.targetId))
			f.write('        <outputPackets>{0}</outputPackets>\n'.format(target.outputPackets))

			if (target.targetId=='SINK'):
				target_map.append(target.targetId)
			else:
				target_map.append(target.targetId)

			target_count.append((target.outputPackets))

			if (target.targetId != 'SINK'):
				gvExpression = "\t%d -> %d [label=\"%d\"]\n"%(int(element.nodeId)+1, int(target.targetId)+1, target.outputPackets)
				gvLines.append(gvExpression)
			for dependency in range(len(target.dataSourcesList)):
				f.write(dependencyBegin)
				if (target.dataSourcesList[dependency]=='LOCAL'):
					f.write('          <sourceType>local</sourceType>\n')
				elif (target.dataSourcesList[dependency]=='SOURCE'):
					f.write('          <sourceType>source</sourceType>\n')
					if element.nodeId not in sourceNodes:
						sourceNodes.append(element.nodeId)
				else:
					f.write('          <sourceType>node</sourceType>\n')
				f.write('          <sourceId>{0}</sourceId>\n'.format(target.dataSourcesList[dependency]))
				f.write('          <triggerPackets>{0}</triggerPackets>\n'.format(target.triggerPacketsList[dependency]))
				f.write(dependencyEnd)

			
				if (target.dataSourcesList[dependency]=='LOCAL' or target.dataSourcesList[dependency]=='SOURCE'):
					if not (target.dataSourcesList[dependency] in src_map):			# Adding SOURCE and LOCAL if not they're not already there
						src_map.append(target.dataSourcesList[dependency])
				else:
					if not (target.dataSourcesList[dependency] in src_map):
						src_map.append(target.dataSourcesList[dependency])

				#if (target.dataSourcesList[dependency]!='LOCAL'):
				#	src_map.append(int(target.dataSourcesList[dependency])-1)
				#else:
				#	src_map.append(element.nodeId)

				#src_count[target_index].append(target.triggerPacketsList[dependency])

			target_index += 1
			#src_count.append([])

			f.write(targetEnd)
		f.write(elementEnd)
		#src_count.pop(-1)

		target_index = 0
		src_map.sort()
		
		for target in element.targetsList:
			srcCts = [-1]*len(src_map)
			for dependency in range(len(target.dataSourcesList)):
				if (target.dataSourcesList[dependency] in src_map):
					srcCts[src_map.index(target.dataSourcesList[dependency])] = target.triggerPacketsList[dependency]
			src_count[target_index] = srcCts
			target_index += 1
			src_count.append([])

		src_count.pop(-1)
					

		# Add the -1s
		while len(target_map) < maxDegree:
			target_map.append('-1')
		while len(target_count) < maxDegree:
			target_count.append(-1)
		while len(src_map) < maxDegree:
			src_map.append('-1')
		while len(src_count) < maxDegree:
			src_count.append([-1]*maxDegree)
		for i in src_count:
			while len(i) < maxDegree:
				i.append(-1)

		# Change the SINK and SOURCE to the actual node IDs

		for ss in range(len(target_map)):
			if target_map[ss] == 'SINK':
				target_map[ss] = element.nodeId

		for ss in range(len(src_map)):
			if src_map[ss] == 'SOURCE':
				src_map[ss] = element.nodeId	
		

		# Writing C code

		c.write("\t//Node %s\n" % (element.nodeId))
		c.write("\t\t//Target Map\n\t\t")
		c.write(", ".join(target_map))
		c.write(',\n')
		c.write("\t\t//Target Count\n\t\t")
		#[ c.write(", ".join(j)) for j in repr(target_count) ]
		c.write(", ".join(map(str, target_count)) )
		c.write(',\n')
		c.write("\t\t//Source Map\n\t\t")
		c.write(", ".join(src_map))
		c.write(',\n')
		c.write("\t\t//Source Count\n\t\t")
		for i in range(len(src_count)-1):
			c.write(", ".join(map(str, src_count[i])) )
			c.write(",\n\t\t")
		c.write(", ".join(map(str, src_count[-1])) )
		#c.write(", ".join(src_count[-1]))
		c.write(",\n")

		for target in element.targetsList:			
			sourceRatios = []
			for srcs in range(len(target.dataSourcesList)):
				coefficient = float(target.outputPackets)/target.triggerPacketsList[srcs]
				for coeffs in networkCoefficients:
					if element.nodeId == coeffs[1] and target.dataSourcesList[srcs] == coeffs[0]:
						sourceRatios.append(coeffs[3]/target.triggerPacketsList[srcs])
			if sourceRatios:
				minIndex = np.argmin(sourceRatios)
				coefficient = sourceRatios[minIndex]*float(target.outputPackets)
			else:
				coefficient = float(target.outputPackets)/target.triggerPacketsList[0] # This takes care of special case of 'SOURCE'
			networkCoefficients.append([element.nodeId, target.targetId, target.outputPackets, coefficient])

	graph.networkBandwidths = networkCoefficients

	c.write("\t-2\n\t};\n")

	#gf.close()

	for j in graph.networkBandwidths:
		print j
	
	f.write(listEnd)
	f.close()


	f = open(filename, 'r')
	data = f.readlines()
	f.close()

	#g.write('\n')
	#for i in data:
	#	g.write(i)
	#g.write('</Graph>\n')
	#g.close()

	f = open(gvName, 'r')
	data = f.readlines();
	f.close()

	editedF = []
	edgeCount = 0
	
	re.DOTALL
	beenDone = False
	

	line = 0
	#while (data[line] != '}'):
	while (not (re.findall(r'^}', data[line]))):
	#for line in range(len(data)):
		trafficLine = re.findall(r'.*->.*', data[line])
		nodeLine = re.findall(r'.*\[label="P.*"', data[line])
		numbers = re.findall(r'\d+', data[line])
		
		if nodeLine and not beenDone:
			l = "  {0} [label=\"P{1}\"]\n".format(numbers[0], repr(int(numbers[1])-1))
			if numbers[1]=='0':
				beenDone = True
				editedF.append(data[line])
			elif int(numbers[0])!=numberNodes:
				editedF.append(l)
			else:
				l = "  {0} [label=\"P{1}\"]\n".format(numbers[0], repr(int(numbers[1])-1))
				l += "  node [margin=0 fontname=arial fontcolor=black fontsize=10 shape=circle width=0.5 fixedsize=true style=filled fillcolor=forestgreen]\n"
				l += "  {0} [label=\"Source\"]\n".format(numberNodes+1)
				l += "  node [margin=0 fontname=arial fontcolor=black fontsize=12 shape=circle width=0.5 fixedsize=true style=filled fillcolor=darksalmon]\n"
				l += "  {0} [label=\"Sink\"]\n".format(numberNodes+2)
				editedF.append(l)
			line += 1
			continue


		elif not trafficLine:
			editedF.append(data[line])
			line += 1
			continue

		elif trafficLine:
			originNode = int(numbers[0])-1
			destinationNode = int(numbers[1])-1
			for el in graph.networkBandwidths:
				if int(el[0]) == originNode and int(el[1]) == destinationNode and el[1]!=el[0]:
					s = '   {0} -> {1} [label="{2:.3}"]\n'.format(numbers[0], numbers[1], el[3])
					editedF.append(s)
					break
					
			line += 1
			continue
		#else:
		#	editedF.append(gvLines[edgeCount])
		#	edgeCount += 1

	editedF.append('}\n')

	for line in range(len(editedF)):
		if re.findall(r'rank=', editedF[line]):
			editedF.insert(line, "  {{rank=same {0}}}\n".format(numberNodes+1))
			editedF.insert(line, '\n')
			for el in graph.networkBandwidths:
				if el[1] == 'SINK':
					s = '   {0} -> {1} [label="{2:.3}"]\n'.format(int(el[0])+1, numberNodes+2, el[3])
					editedF.insert(line, s)
			for src in sourceNodes:
				s = '   {0} -> {1} [label="{2:.3}"]\n'.format(numberNodes+1, int(src)+1, 1.0)
				editedF.insert(line,s)
			break

	for line in range(len(editedF)-1, 0, -1):
		if re.findall(r'rank=', editedF[line]):
			editedF.insert(line+1, "  {{rank=same {0}}}\n".format(numberNodes+2))
			break

	f = open(gvName, 'w')
	for line in editedF:
		f.write(line)
	f.close() 

	os.system("dot %s -Tpdf -o %s"%(gvName, gvName[:-3]+'.pdf'))
	
	

def fileSeek(directory, extension='xml', extensionLength=3):
	fileList = []		
	for f in os.listdir(directory):
		if f[-extensionLength:] == extension:
			fileList.append(f)
	return fileList


def main(*args):
	try:

		#fileList = fileSeek('Archive/', 'xml', 3)
		seed = random.seed(100)
		np.random.seed(100)

		graphID = 0
		fileList, fL = [], os.listdir("Archive")
		for i in fL:
			if i[-4:] == ".xml":
				fileList.append(i)

		for i in fileList:

			fOutName = 'Traffic/'+i
			fInName  = 'Archive/'+i
			gvName	 = 'Traffic/'+i[:-4] +'.gv'
			graphID  = int(i[3:-4])

			# Import XML structure into Python classes
			#graph = importXML(fInName, graphID)
			graph, dumpData = importXML(fInName, graphID)
			# Create random traffic based on imported graphs
			print "Graph ID:", graphID
			trafficGraph = addRandomTraffic(graph, graphID, connectivity=0.7, dependenceConstant=0, localThreshold=0.8, inputPacketRange=8, outputPacketRange=5, maxDegree=8)
			# Create XML for traffic pattern
			exportXML(fOutName, trafficGraph, dumpData, gvName)
			
		# For every graph in folder
		#for filename in range(len(fileList)):
		#	fname = 'Archive/'+fileList[filename]
			# Import XML structure into Python classes
		#	graph, dumpData = importXML(fname, 1)
			# Create random traffic based on imported graphs
		#	trafficGraph = addRandomTraffic(graph, filename+1, connectivity=0.7, dependenceConstant=10, localThreshold=0.8, inputPacketRange=8, outputPacketRange=5)
			# Create XML for traffic pattern
		#	exportXML(fileList[filename], trafficGraph, dumpData, gvName)
		

	except:
		print "Unexpected error:", sys.exc_info()[0]
		raise

	else:
		return 0

if __name__ == '__main__':
	sys.exit(main(*sys.argv))
