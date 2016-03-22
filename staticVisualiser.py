################################# Traffic Generator v1.0 ################################
#											#
# 				    University of York					#
# 				     Graceful Project					#
# 				  Created by Pedro Campos				#
# 					York, 2015					#
#											#
#########################################################################################

# -*- coding: utf-8 -*-

import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import re
import os, errno, sys
import random
import math
import string
import cPickle as pickle
import time
from random import randint, choice
from matplotlib.collections import PatchCollection
import imp
from matplotlib.path import Path
import matplotlib.image as mpimg
from PIL import Image
from operator import itemgetter

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

from graphGeneratorXML import Node
from graphGeneratorXML import Edge
from graphGeneratorXML import TrafficEdge
from graphGeneratorXML import TrafficTarget
from graphGeneratorXML import TrafficElement
from graphGeneratorXML import TrafficList
from graphGeneratorXML import basicGraph
from graphGeneratorXML import Graph

import images2gif

cellSize = 0.6
numberTasks = 9
numberNodes = 9
gridSize = 3
rowSize = 4
colSize = 4
arrowSeparation = 0.1
refreshTime = 0.5
maximumMode = 'relative'
bufWidth = 0.08
bufHeight = 0.025
localWidth = 0.12
localHeight = 0.06
maxIdleTime = 200
linkMax = 10E6
maxPow = 300 # in milliwatts
maxTemp = 90
graphName = 'dag0020'

addRouting = True

degreeSymbol = unichr(176)

codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         ]

def handle_close(evt):
    print('Closed Figure!')

def label(xy, text):
	global cellSize
	y = xy[1] + cellSize/2-0.002  # shift y-value for label so that it's centered with cell
	x = xy[0] + cellSize/2-0.002  # shift x-value for label so that it's centered with cell
	plt.text(x, y, text, ha="center", family='sans-serif', size=18)

def interpretBufferLevels(word):

	coded = []
	interpreted = []
	for i in word:
		#print i, int(i, 16)
		a = bin(int(i, 16) & 0x7)[2:].zfill(4) 
		cleanStr = [False for k in range(len(a)-1)]
		for i in range(1, len(a), 1):
			if a[i] == '0':
				cleanStr[i-1] = False
			else:
				cleanStr[i-1] = True

		if cleanStr[0] and cleanStr[1] and cleanStr[2]:				# All flags active - setup stage
			bufColour = 'black'
		elif (not cleanStr[0] and not cleanStr[1] and not cleanStr[2]):		# No flags active - normal operation
			bufColour = 'green'
		elif cleanStr[2]:							# Empty flag active - empty FIFO
			bufColour = 'blue'
		elif cleanStr[0] or cleanStr[1]:
			bufColour = 'red'
		interpreted.append(bufColour)
		#print cleanStr
	
	return interpreted

def onpick1(event):
	print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            event.button, event.x, event.y, event.xdata, event.ydata)


def main(*args):
	try:
		global cellSize, bufHeight, bufWidth, localWidth, localHeight
		global maximumMode
		global codes
		global maxIdleTime
		global adaptiveRouting
		random.seed(100)

		if not os.path.exists('frames'):
			os.mkdir('frames')

		if 'saveframe' in args:
			saveFrames = True
		else:
			saveFrames = False

		if os.path.exists('frames') and saveFrames:
			os.system('rm ./frames/*.png')

		if 'absolute' in args:
			maximumMode = 'absolute'
			print "Enter maximum throughput for a link:"
			maxThroughput = float(input())
			maxNodeThroughput = 5*maxThroughput
		else:
			maximumMode = 'relative'

		if 'balanced' in args:
			adaptiveRouting = True
		else:
			adaptiveRouting = False

		#######################################################################
		#		Adjustments for multi-grid plotting
		#######################################################################

		cellSize = cellSize/gridSize
		#arrowSeparation = 0.1*cellSize

		#######################################################################
		#		Adjustments for multi-grid plotting
		#######################################################################


		# Interpreting graph from data file

		if numberNodes == 100:
			f = open("../100nodes/Traffic/{0}.data".format(graphName), "r")
		elif numberNodes == 1000:
			f = open("../1000nodes/Traffic/{0}.data".format(graphName), "r")
		else:
			f = open("../Traffic/{0}.data".format(graphName), "r")
		#f = open("Traffic/dag0001.data", "r")
		graph = pickle.load(f)
		f.close()			

		# Set up figure and axes
		
		#fig, ax = plt.subplots(figsize=(40,24))
		fig, ax = plt.subplots(figsize=(24,24))
		ax.set_picker(True)
		fig.set_picker(True)
		fig.canvas.mpl_connect('close_event', handle_close)
		#fig = plt.figure(figsize=(20,12))
		#ax = plt.subplots()
		fig.canvas.mpl_connect('pick_event', onpick1)

		#fig.suptitle('Graceful Traffic Visualiser', fontsize=18, x=0.78, y=0.08, ha='right', family='sans-serif')
		plt.axis('off')

		plt.ion()
		plt.show()

		start = time.time()

		fname = "taskMapping.txt"
		ppp = 'touch ' + fname
		os.system(ppp)

		# From this point on, files are scanned and the graph is updated.

		iteration = 1

		stopPlotting = False
		activeFile = True

		while(1):
			ax.set_picker(True)
			fig.set_picker(True)

			# If nothing happens for 10 seconds, leave
			if (time.time() - os.path.getmtime(fname) > maxIdleTime):	
				#stopPlotting = False
				activeFile = False
			else:
				activeFile = True
			stopPlotting = not activeFile

			# Interpreting task mapping from text file, just to know where to write tasks

			'''taskMapping = np.ones(numberNodes, dtype=int)
			taskMapping = -taskMapping

			mapData = open('taskMapping.txt','r').readlines()
			for line in range(len(mapData)):
				#mapLine = re.findall(r"\$\$ .*", mapData[line])
				headerLine = re.findall(r"NODE", mapData[line])
				numbers = re.findall(r"[+-]?\d+", mapData[line])
				endLine = re.findall("{\#}*", mapData[line])
				if numbers and not headerLine:
					taskMapping[int(numbers[1])] = int(numbers[0])
				
				if endLine:
					break'''

			#taskMapping = [0, 1, 2, 23, 4, 15, 6, 31, 8]
			#taskMapping = [-1, -1, -1, -1, -1, -1, -1, -1, -1]
		
			#taskMapping = [0, 1, 2, 3, 4, -1, -1, -1, -1]

			taskMapping = [[566514001, 0], [594308093, 1], [979738199, 2], [818607243, 3], [142299675, 4], [435137665, 5], [583129108, 6], [357174482, 7], [45719769, 8], [583129168, 9], [357184482, 10], [457769, 11]]

			taskMapping.sort(key=itemgetter(0), reverse=int(0.0))

			taskMapping = taskMapping[0:numberTasks]
			taskMapping = [element[1] for element in taskMapping]
			print taskMapping


			#########################################################

			# Overlaying deterministic routing traffic on figure

			#########################################################

			networkCoefficients = []

			if(addRouting):

				connectMatrix = np.ndarray(shape=(len(graph.trafficElements),len(graph.trafficElements)), dtype=int)

				connectMatrix[:,:] = 0

				sinkedElements = []

				for element in graph.trafficElements:
					#print "In node ", element.nodeId, ":"
					for target in element.targetsList:
						if target.targetId != 'SINK':
							connectMatrix[int(element.nodeId), int(target.targetId)] = target.outputPackets
						else:
							sinkedElements.append([int(element.nodeId), target.outputPackets])

						sourceRatios = []

						for srcs in range(len(target.dataSourcesList)):
							coefficient = float(target.outputPackets)/target.triggerPacketsList[srcs]
							for coeffs in networkCoefficients:
								if element.nodeId == coeffs[1] and target.dataSourcesList[srcs] == coeffs[0]:
									sourceRatios.append(coeffs[3]/target.triggerPacketsList[srcs])
									

						if sourceRatios:
							minIndex = np.argmin(sourceRatios)
							coefficient = sourceRatios[minIndex]*float(target.outputPackets)
							#print "Node ID: ", element.nodeId, "Target ID: ", target.targetId, "Output packets:", target.outputPackets, "Min Index: ", minIndex, "Ratio: ", float(target.outputPackets)/target.triggerPacketsList[minIndex], "New ratio: ", coefficient
							
						else:
							coefficient = float(target.outputPackets)/target.triggerPacketsList[0] # This takes care of special case of 'SOURCE'

						networkCoefficients.append([element.nodeId, target.targetId, target.outputPackets, coefficient])		

						#for srcs in sourceRatios:
						#	networkCoefficients.append(srcs)

				print "**************"

				print "Network Coefficients:"
				for el in networkCoefficients:
					print el

				print "&&&&&&&&&&&&&&"

				# Create routing matrix

				routingMatrix = np.zeros((len(graph.trafficElements), len(graph.trafficElements)), dtype=int)

				edgeList = []
	
				for i in range(len(connectMatrix)):
					for j in range(len(connectMatrix[i])):
						if connectMatrix[i,j] != 0:
							edgeList.append([taskMapping[i], taskMapping[j], connectMatrix[i,j]])

				hopList = []

				localTraffic2PS = np.zeros((rowSize, colSize))
				graphEnergy = 0

				for i in range(len(edgeList)):

					originNode = edgeList[i][0]
					destinationNode = edgeList[i][1]
					packets = edgeList[i][2]

					destCol = destinationNode%colSize
					originCol = originNode%colSize
					destRow = int(float(destinationNode)/colSize)
					originRow = int(float(originNode)/colSize)
					print originRow, originCol, destRow, destCol
			
					if originNode != -1 and destinationNode != -1:
						hopList.append([originRow, originCol, destRow, destCol, destCol-originCol, destRow-originRow, packets, originNode, destinationNode, networkCoefficients[i][3]])
						graphEnergy += networkCoefficients[i][3]
						localTraffic2PS[destRow, destCol] += networkCoefficients[i][3]


				# Iterate through connections and update link load

				maxCol = colSize-1
				maxRow = rowSize-1

				posHorizontalTraffic = np.zeros((rowSize, maxCol))
				posVerticalTraffic = np.zeros((maxRow, colSize))
				negHorizontalTraffic = np.zeros((rowSize, maxCol))
				negVerticalTraffic = np.zeros((maxRow, colSize))

				#########################################################

				# 	Performing deterministic x-y routing

				#########################################################
			
				pathRigidness = 0
				totalBandwidth   = 0
				maxMultiplicity = (math.factorial(maxCol+maxRow))/(math.factorial(maxCol)*math.factorial(maxRow))

				if not adaptiveRouting:
					for i in hopList:
						print "hop", i

						# Route through x-orientation
						originRow = i[0]
						destRow = i[2]
						originCol = i[1]
						destCol = i[3]
						colDiff = i[4]
						rowDiff = i[5]
						packets = i[9]

						pathMultiplicity = (math.factorial(abs(colDiff)+abs(rowDiff)))/(math.factorial(abs(colDiff))*math.factorial(abs(rowDiff)))
						pathRigidness += (float(pathMultiplicity)/maxMultiplicity * packets)
						totalBandwidth += packets

						if (colDiff != 0):
							if(colDiff > 0):
								for j in range(colDiff):
									posHorizontalTraffic[originRow, originCol+j] += packets
							else:
								for j in range(abs(colDiff)):
									negHorizontalTraffic[originRow, originCol-1-j] += packets

	
						# Route through y-orientation
						if (rowDiff != 0):
							if(rowDiff > 0):
								for j in range(rowDiff):
									posVerticalTraffic[originRow+j, destCol] += packets # Destination column used because x-routing has already been performed
							else:
								for j in range(abs(rowDiff)):
									negVerticalTraffic[originRow-1-j, destCol] += packets # Destination column used because x-routing has already been performed

				#########################################################

				# 	Path-based routing, extracting information

				#########################################################
				
				elif adaptiveRouting:

					paths = []
	
					random.shuffle(hopList)
			
					for i in hopList:
						print "hop", i
					
						# Route through x-orientation
						originRow = i[0]
						destRow = i[2]
						originCol = i[1]
						destCol = i[3]
						colDiff = i[4]
						rowDiff = i[5]
						packets = i[9]

						currRow, currCol = originRow, originCol

						pathMultiplicity = (math.factorial(abs(colDiff)+abs(rowDiff)))/(math.factorial(abs(colDiff))*math.factorial(abs(rowDiff)))
						pathRigidness += (float(pathMultiplicity)/maxMultiplicity * packets)
						totalBandwidth += packets

						while colDiff!=0 or rowDiff!=0:

							print "Path multiplicity:", pathRigidness
					
							if colDiff > 0:
								if rowDiff > 0:
									if posHorizontalTraffic[currRow,currCol] > posVerticalTraffic[currRow,currCol]:
										posVerticalTraffic[currRow,currCol] += packets
										currRow += 1
									else:
										posHorizontalTraffic[currRow,currCol] += packets
										currCol += 1
								elif rowDiff < 0:
									if posHorizontalTraffic[currRow,currCol] > negVerticalTraffic[currRow-1,currCol]:
										negVerticalTraffic[currRow-1,currCol] += packets
										currRow -= 1
									else:
										posHorizontalTraffic[currRow,currCol] += packets
										currCol += 1
								else:
									posHorizontalTraffic[currRow,currCol] += packets
									currCol += 1
								
							elif colDiff < 0:
								if rowDiff > 0:
									if negHorizontalTraffic[currRow,currCol-1] > posVerticalTraffic[currRow,currCol]:
										posVerticalTraffic[currRow,currCol] += packets
										currRow += 1
									else:
										negHorizontalTraffic[currRow,currCol-1] += packets
										currCol -= 1
								elif rowDiff < 0:
									if negHorizontalTraffic[currRow,currCol-1] > negVerticalTraffic[currRow,currCol-1]:
										negVerticalTraffic[currRow,currCol-1] += packets
										currRow -= 1
									else:
										negHorizontalTraffic[currRow,currCol-1] += packets
										currCol -= 1
								else:
									negHorizontalTraffic[currRow,currCol-1] += packets
									currCol -= 1

							else:
								if rowDiff > 0:
									posVerticalTraffic[currRow,currCol] += packets
									currRow += 1
								elif rowDiff < 0:
									negVerticalTraffic[currRow,currCol-1] += packets
									currRow += 1
							

							colDiff = destCol-currCol
							rowDiff = destRow-currRow

				print "Total bandwidth:", totalBandwidth

				pathRigidness = 2*(totalBandwidth/(totalBandwidth+pathRigidness))-1
				#pathMultiplicities.append(pathRigidness)

				print "Average path multiplicity:", pathRigidness	

				maxPosX = np.amax(posHorizontalTraffic.reshape(1,-1))
				maxNegX = np.amax(negHorizontalTraffic.reshape(1,-1))
				maxPosY = np.amax(posVerticalTraffic.reshape(1,-1))
				maxNegY = np.amax(negVerticalTraffic.reshape(1,-1))

				#posMaxTraffic = np.amax(np.maximum(posVerticalTraffic.reshape(1,-1), posHorizontalTraffic.reshape(1,-1)))
				#negMaxTraffic = np.amax(np.maximum(negVerticalTraffic.reshape(1,-1), negHorizontalTraffic.reshape(1,-1)))

				maxLinkTraffic = np.amax([maxPosX, maxPosY, maxNegX, maxNegY])

				if maximumMode == 'relative':
					maxThroughput = maxLinkTraffic

				# Create nxn mesh to store traffic through each node
				nodeTraffic = np.zeros(shape=(rowSize,colSize))

				networkEnergy = 0

				for i in range(maxCol):
					for j in range(rowSize):

						networkEnergy += posHorizontalTraffic[j,i]
						networkEnergy += negHorizontalTraffic[j,i]

						nodeTraffic[j, i] += posHorizontalTraffic[j,i]
						nodeTraffic[j, i+1] += 2*posHorizontalTraffic[j,i]
						nodeTraffic[j, i+1] += negHorizontalTraffic[j,i]
						nodeTraffic[j, i] += 2*negHorizontalTraffic[j,i]

				for i in range(maxRow):
					for j in range(colSize):

						networkEnergy += posVerticalTraffic[i,j]
						networkEnergy += negVerticalTraffic[i,j]
				
						# Add traffic to both the outgoing node and the incoming one
						nodeTraffic[i, j] += posVerticalTraffic[i,j]
						nodeTraffic[i+1, j] += 2*posVerticalTraffic[i,j]
						nodeTraffic[i+1, j] += negVerticalTraffic[i,j]
						nodeTraffic[i, j] += 2*negVerticalTraffic[i,j]
	
				print "#################"

				print posHorizontalTraffic
				print negHorizontalTraffic
				print posVerticalTraffic
				print negVerticalTraffic

			#########################################################

			# Finished overlaying deterministic routing traffic on figure

			#########################################################
			
			# create nxn grid to plot the artists
			grid = np.mgrid[1.05:0.05:rowSize*1j, 0.05:1.05:colSize*1j].reshape(2, -1).T

			# Adjusting origin of graph to be at the top left corner
			grid[:, 0], grid[:, 1] = grid[:, 1], grid[:, 0].copy()

			rectSize = 1.05 + cellSize

			backXPos = grid[-rowSize][0] - cellSize/2
			backYPos = grid[-colSize][1] - cellSize/2

			backXSize = grid[-1][0] - grid[rowSize][0] + 2*cellSize
			backYSize = grid[0][1] - grid[-1][1] + 2*cellSize

			rect = mpatches.Rectangle([backXPos, backYPos], backXSize, backYSize, ec="none", color='grey', alpha=0.1)
			#ax.add_patch(rect)

			#ax.text(grid[j*gridSize+i, 0]+cellSize, grid[j*gridSize+i, 1]+(0.5+arrowSeparation)*cellSize, '{0}'.format(int(negHorizontalTraffic[j,i])), fontsize=14, ha='center', family='sans-serif')
			#ax.text(grid[j*gridSize+i, 0]+cellSize*2.5, grid[j*gridSize+i, 1]+(0.5-arrowSeparation)*cellSize, '{0}'.format(int(posHorizontalTraffic[j,i])), fontsize=14, ha='center', family='sans-serif')
			#ax.text(grid[i*gridSize+j, 0]+(0.5+arrowSeparation)*cellSize, grid[i*gridSize+j, 1], '{0}'.format(int(posVerticalTraffic[i,j])), fontsize=14, ha='center', family='sans-serif')
			#ax.text(grid[i*gridSize+j, 0]+(0.5-arrowSeparation)*cellSize, grid[i*gridSize+j, 1]-1.5*cellSize, '{0}'.format(int(negVerticalTraffic[i,j])), fontsize=14, ha='center', family='sans-serif')

			
			maxWidth = 0.1
			arrowLength = grid[1][0] - grid[0][0] - cellSize

			arrows = []
			arrowColours = []

			#  Drawing horizontal edges between nodes

			for i in range(maxCol):
				for j in range(rowSize):

					# Draw horizontal lines

					if posHorizontalTraffic[j,i] != 0:
						fillBool = True
					else:
						fillBool = False

					arrow = mpatches.Arrow(grid[j*colSize+i, 0]+cellSize, grid[j*colSize+i, 1]+(0.5+arrowSeparation)*cellSize, arrowLength, 0, width=0.2/max(colSize, rowSize), fill=fillBool, color='white', visible=fillBool)

					if fillBool:
						arrows.append(arrow)
						arrowColours.append(posHorizontalTraffic[j,i]/maxThroughput)

					if negHorizontalTraffic[j,i] != 0:
						fillBool = True
					else:
						fillBool = False

					arrow = mpatches.Arrow(grid[j*colSize+i, 0]+cellSize+arrowLength, grid[j*colSize+i, 1]+(0.5-arrowSeparation)*cellSize, -arrowLength, 0, width=0.2/max(colSize, rowSize), color='white', fill=fillBool, visible=fillBool)
					if fillBool:
						arrows.append(arrow)
						arrowColours.append(negHorizontalTraffic[j,i]/maxThroughput)

			# Drawing vertical edges between nodes

			for i in range(maxRow):
				for j in range(colSize):
					
					# Draw vertical arrows

					if posVerticalTraffic[i,j] != 0:
						fillBool = True
					else:
						fillBool = False

					arrow = mpatches.Arrow(grid[i*colSize+j, 0]+(0.5+arrowSeparation)*cellSize, grid[i*colSize+j, 1], 0, -arrowLength, width=0.2/max(colSize, rowSize), fill=fillBool, color='white', visible=fillBool)
					if fillBool:
						arrows.append(arrow)
						arrowColours.append(posVerticalTraffic[i,j]/maxThroughput)

					if negVerticalTraffic[i,j] != 0:
						fillBool = True
					else:
						fillBool = False
					arrow = mpatches.Arrow(grid[i*colSize+j, 0]+(0.5-arrowSeparation)*cellSize, grid[i*colSize+j+colSize, 1]+cellSize, 0, arrowLength, width=0.2/max(colSize, rowSize), fill=fillBool, color='white', visible=fillBool)
					if fillBool:
						arrows.append(arrow)
						arrowColours.append(negVerticalTraffic[i,j]/maxThroughput)


			#print "local traffic", localTraffic2PS[i,j]
			#nodeTraffic = [[(nodeTraffic[i,j] + localTraffic2PS[i,j]) for j in range(len(nodeTraffic[i]))] for i in range(len(nodeTraffic))]
			#print "node traffic", nodeTraffic[i,j]	

			for i in range(len(nodeTraffic)):
				for j in range(len(nodeTraffic[i])):
					nodeTraffic[i,j] += localTraffic2PS[i,j]			

			nodes = []
			nodeColours = []
			emptyNodes = []

			maxNodeTraffic = np.amax(nodeTraffic)

			if maximumMode == 'relative':
				maxNodeThroughput = maxNodeTraffic

			try: 
				minNodeTraffic = np.amin(nodeTraffic[np.nonzero(nodeTraffic)])
			except (Exception):
				print "No traffic."
				minNodeTraffic = 1

			for i in range(rowSize):
				for j in range(colSize):

					if nodeTraffic[i,j] == 0:
						fillBool = False
						fc = None
					else:
						fillBool = True
						fc = 'white'

					rect = mpatches.Rectangle(grid[i*colSize+j], cellSize, cellSize, fill=fillBool, ec=None, visible=False, color=fc)

					if fillBool:
						nodes.append(rect)
						nodeColours.append(nodeTraffic[i,j]/maxNodeThroughput)
					else:
						emptyNodes.append(rect)

					nodeID = i*colSize+j
					for k in range(len(taskMapping)):
						if taskMapping[k]==nodeID:
							#label(grid[i*gridSize+j], "P{0}".format(k))

							# Write task number inside box
							x = grid[i*colSize+j][0] + cellSize/2
							y = grid[i*colSize+j][1] + cellSize/2

							ax.text(x,y, 'P{0}'.format(k), fontsize=30/pow(max(rowSize,colSize),0.5), ha='center', va='center', family='sans-serif')

							#circle = mpatches.Circle([x,y], radius=cellSize/4, fill=False)
							#ax.add_patch(circle)

							#plt.savefig('testRoute.png',bbox_inches='tight')

			if arrows:
				collection = PatchCollection(arrows, cmap='Oranges', alpha=0.99)
				collection.set_array(np.array(arrowColours))
				ax.add_collection(collection)
				#norm = matplotlib.colors.Normalize(vmin=minNodeTraffic, vmax=maxNodeTraffic)
				#plt.colorbar(collection, shrink=0.8, norm=norm)
				#pcm = ax.pcolor(grid, norm=matplotlib.colors.Normalize(vmin=minNodeTraffic, vmax=maxNodeTraffic), cmap='Oranges', alpha=0.8)
				#fig.colorbar(pcm, ax=ax, shrink=0.2)

			if nodes:
				collection = PatchCollection(nodes, cmap='Oranges', alpha=0.99, match_original=False)
				collection.set_array(np.array(nodeColours))
				ax.add_collection(collection)
				#plt.colorbar(collection)
				#collection.set_visible(False)

			if emptyNodes:
				collection = PatchCollection(emptyNodes, alpha=0.8, match_original=True)
				ax.add_collection(collection)

			#rect = mpatches.Rectangle([backXPos, backYPos], backXSize, backYSize, ec="none", color='grey', alpha=0.1)		
			


			#cax = fig.add_axes([1.40, 0.05, 0.035, 0.09], label='axes1', autoscale_on=True)
			#norm = matplotlib.colors.Normalize(vmin=minNodeTraffic, vmax=maxNodeTraffic)

			#cb1 = matplotlib.colorbar.ColorbarBase(cax, norm=norm, orientation='vertical', cmap='Oranges')
				


			# Final commands

			plt.axis('equal')
			plt.axis('off')
			filename = "frames/staticFrame" + repr(iteration) + ".png"
			plt.draw()
			if saveFrames:		
				plt.savefig(filename,bbox_inches='tight')
			time.sleep(refreshTime)
			plt.cla()
			plt.clf()
			ax = fig.add_subplot(111)
			iteration += 1

			if stopPlotting:
				print "A total of", maxIdleTime, "seconds went by without activity in any file."
				break
			

	except:
		print "Unexpected error:", sys.exc_info()[0]
		raise

	else:
		return 0

if __name__ == '__main__':
	sys.exit(main(*sys.argv))
