#!/usr/bin/env python
#  analysisTools.py
#  Created by nicain on 11/4/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.

################################################################################
# This function plots creates a multiline speed accuracy plot:
def speedAccuracyMultiLine(sliceDict, saveResultDir = './savedResults', whichRun = 0, tDel = 2000, tPen = 0, tND = 350, newFigure = 1, quickName = -1, N = 5, colorBar = 1, titleString = -1, lims = -1, plotLabels = 1, color = [], saveFig=0):
	from numpy import array, linspace, inf
	from pylab import figure, subplots_adjust, cm, flipud, pcolor, colorbar, hold, savefig
	import pylab as pl
	import copy
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	seqDimensionTuple = sliceDict['MultiLineVar']
	if isinstance(seqDimensionTuple, str):
		seqDimension = seqDimensionTuple
		settings, FD, numberOfJobs, gitVersion =  getSettings(quickName, saveResultDir, whichRun=whichRun)
		vals = settings[seqDimension]
		seqDimensionList = linspace(min(vals), max(vals), N)
	else:
		seqDimension = seqDimensionTuple[0]
		seqDimensionList = seqDimensionTuple[1]
		
	if color == []:
		seqDimensionListArray = array(seqDimensionList, dtype=float)
		colorMatrix=cm.autumn_r((seqDimensionListArray-min(seqDimensionListArray))/(max(seqDimensionListArray)-min(seqDimensionListArray)))
	if colorBar: pcolor(array([[min(seqDimensionList),max(seqDimensionList)]]),cmap=cm.autumn_r,visible=False)

	del sliceDict['MultiLineVar']
	if newFigure:
		figure()
	minX = inf
	minY = inf
	maxX = -inf
	maxY = -inf
	figure(99)
	for i in range(len(seqDimensionList)):
		sliceDict[seqDimension] = seqDimensionList[i]
		thisPlot = speedAccuracy(copy.copy(sliceDict),saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, newFigure = 0, plotLabels = -1, lims = -1, color = colorMatrix[i])
		if min(thisPlot[0].get_ydata()) < minY:
			minY = min(thisPlot[0].get_ydata())
		if max(thisPlot[0].get_ydata()) > maxY:
			maxY = max(thisPlot[0].get_ydata())
		if min(thisPlot[0].get_xdata()) < minX:
			minX = min(thisPlot[0].get_xdata())
		if max(thisPlot[0].get_xdata()) > maxX:
			maxX = max(thisPlot[0].get_xdata())
	pl.close(99)
			
	if lims == -1:
		lims = ((minX, maxX),(minY, maxY))
	
	if newFigure:
		pl.figure()
	else:
		figure(1)		
		
	for i in range(len(seqDimensionList)):
		sliceDict[seqDimension] = seqDimensionList[i]
		if titleString == -1:
			titleString = ''
		thisPlot = speedAccuracy(copy.copy(sliceDict),saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, lims = lims, newFigure = 0, plotLabels = -1, color = colorMatrix[i])	
		
	if colorBar: 
		cb = colorbar()
		cb.set_label('Color variable: ' + seqDimension)
		
	if plotLabels == 1:
		pl.xlabel('RT')
		pl.ylabel('FC')
		
	if saveFig != 0:
		savefig('/Users/Nick/Desktop/fig1.eps')
	
	

	return lims

################################################################################
# This function plots speed accuracy tradeoff function:
def speedAccuracy(sliceDict, saveResultDir = './savedResults',tDel = 2000, tPen = 0, tND = 350, whichRun = 0, newFigure = 1, quickName = -1, saveFig=0, plotLabels = 1, color = -1, lims = -1,pade=0, makePlot=1, thetaN = 'Default'):
	import pylab as pl
	from numpy import transpose, shape, squeeze, ndarray, array, mean, exp
	import copy

	# Import necessary stuff for pade approx:
	import scipy as sp
	from scipy.optimize import leastsq
	from random import random as r
	
	# Define stuff for Pade:
	maxRecursion = 15
	tol=.15

	#---------------------------------------------------------------------------
	def padeApproxOfY(x,y, N='Default', level=maxRecursion):
		
		# Escape condition:
		if level == 0:
			print 'Max Level Reached.  Utter Failure!'
			return 1,1
			
		# Interpret length of x-vals being returned:
		if N == 'Default':
			N=len(x)
		
		# Set up functions:
		def res(p,y,x): return y-peval(x,p)
		def peval(x,p): return ((p[0]+p[1]*x+p[2]*x**2+p[3]*x**3+p[4]*x**4)/
									(p[5]+p[6]*x+p[7]*x**2+p[8]*x**3+p[9]*x**4))
		
		# Optimize:
		if level == maxRecursion:
			p0=[1,1,1,1,1,1,1,1,1,1]
		elif level == maxRecursion-1:
			p0=[317.41537167, -1146.79958272, 1695.42313448, 715.40164079, 79.73275651,
			    462.86872028, -1719.11197335, 2601.63857965, 839.88221669, 90.31055117]
		else:
			p0 = [r(),r(),r(),r(),r(),r(),r(),r(),r(),r()]
		plsq = leastsq(res, p0, args=(y,transpose(x)),maxfev=5000)
		print plsq[0]
		
		# Return interpolated values:
		def f(x): return peval(array(x),plsq[0])
		newX = sp.linspace(min(x),max(x),N)
		smoothY = f(newX)
		
		if (plsq[1] != 1) or (pl.linalg.norm(y/max(y)-f(x)/max(f(x))) > tol):
			newLevel = level-1
			newX, smoothY = padeApproxOfY(x,y, N=N, level=newLevel)
		
		return newX, smoothY
	#---------------------------------------------------------------------------

	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')
	
	# Get data:
	crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion =  getDataAndSettings(quickName=quickName, saveResultDir=saveResultDir, whichRun=whichRun)
	crossTimeData += tND 
    
	# Reorder dimension list and cube to put theta variable last:
	if 'theta' in sliceDict:
		permuteList = range(len(dims))
		whereIsXDim = dims.index('theta')
		dims[-1], dims[whereIsXDim] = dims[whereIsXDim], dims[-1]
		permuteList[-1], permuteList[whereIsXDim] = permuteList[whereIsXDim], permuteList[-1]
		crossTimeData = transpose(crossTimeData,permuteList)
		resultData = transpose(resultData,permuteList)
	
	# Collapse all non-constant dimensions:
	crossDims = dims[:]
	resultDims = dims[:]
	for collapseDim in iter(sliceDict):
		crossTimeData, resultData, crossDims = reduce1D(crossTimeData, resultData, crossDims, collapseDim, settings[collapseDim], sliceDict[collapseDim], tDel = tDel, tPen = tPen, tND = tND)
	crossTimeSlice = squeeze(crossTimeData)
	resultSlice = squeeze(resultData)
	
	# Get x and y data:
	sliceDict['XVar'] = 'theta'
	counter = -1
	X=[0]*2
	Y=[0]*2
	for yAxis in ['RT','FC']:
		counter += 1
		# Get data:
		X[counter], Y[counter] = (export1D(copy.copy(sliceDict), yAxis, whichRun=whichRun, 
					  quickName=quickName, saveResultDir=saveResultDir))
	
	if pade:
		if thetaN == 'Default':
			thetaN=len(X[0])
		thetaVals, myRT = padeApproxOfY(X[0],Y[0], N=thetaN)
		pl.figure(201)
		pl.plot(thetaVals, myRT)
		pl.plot(X[0],Y[0])

		thetaVals, myFC = padeApproxOfY(X[1],Y[1], N=thetaN)
		pl.figure(202)
		pl.plot(thetaVals, myFC)
		pl.plot(X[1],Y[1])
	else:
		myRT, myFC = Y[0],Y[1]
		thetaVals = X[0]


	
	if makePlot == 1:
		if newFigure:
			pl.figure()
		
		if not(isinstance(color,ndarray)):
			myPlot = pl.plot(myRT,myFC)
			if pade:
				pl.plot(Y[0],Y[1])
		else:
			myPlot = pl.plot(myRT,myFC, color = color)

		if lims != -1:
			pl.xlim(lims[0][0], lims[0][1])
			pl.ylim(lims[1][0], lims[1][1])
			
		if plotLabels == 1:
			pl.xlabel('RT')
			pl.ylabel('FC')
		
		if saveFig != 0:
			pl.savefig('/Users/Nick/Desktop/fig1.eps')
	
	return thetaVals,myRT,myFC



################################################################################
# This function plots a sequence  of 1-D multi plots:
def plot1DSeqMultiLine(sliceDict, whatToPlot, saveResultDir = './savedResults', whichRun = 0, tDel = 2000, tPen = 0, tND = 350, newFigure = 1, quickName = -1, seqLength = 4, N=5, saveFig=0, colorBar=1):
	from numpy import array, linspace, inf
	from pylab import figure, subplot, suptitle, subplots_adjust, savefig, ylim
	import copy
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')
	if whatToPlot == 'All':
		plot1DSeqMultiLine(copy.copy(sliceDict), 'RR', saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, newFigure = newFigure, quickName = quickName, seqLength = seqLength, N = N)
		plot1DSeqMultiLine(copy.copy(sliceDict), 'RT', saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, newFigure = newFigure, quickName = quickName, seqLength = seqLength, N = N)
		plot1DSeqMultiLine(copy.copy(sliceDict), 'FC', saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, newFigure = newFigure, quickName = quickName, seqLength = seqLength, N = N)
		return

	seqDimensionTuple = sliceDict['SeqVar']
	if isinstance(seqDimensionTuple, str):
		seqDimension = seqDimensionTuple
		settings, FD, numberOfJobs, gitVersion =  getSettings(quickName, saveResultDir, whichRun = whichRun)
		vals = settings[seqDimension]
		seqDimensionList = linspace(min(vals), max(vals), seqLength)
	else:
		seqDimension = seqDimensionTuple[0]
		seqDimensionList = seqDimensionTuple[1]

	del sliceDict['SeqVar']
	if newFigure:
		figure(num=None,figsize=(4*seqLength, 4))
	minY = inf
	maxY = -inf
	for i in range(len(seqDimensionList)):
		subplot(1,len(seqDimensionList),i+1)
		sliceDict[seqDimension] = seqDimensionList[i]
		titleString = seqDimension + '=' + '%-5.3f' % seqDimensionList[i]
		yLimsBack = plot1DMultiLine(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, N = N, newFigure = 0, plotYLabel = 0, yLims = -1, colorBar = 0)
		if yLimsBack[0] < minY:
			minY = yLimsBack[0]
		if yLimsBack[1] > maxY:
			maxY = yLimsBack[1]
			
	yLims = (minY, maxY)
			
	for i in range(len(seqDimensionList)):
		subplot(1,len(seqDimensionList),i+1)
		sliceDict[seqDimension] = seqDimensionList[i]
		titleString = seqDimension + '=' + '%-5.3f' % seqDimensionList[i]
		if i==0:
			thisPlot = plot1DMultiLine(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, yLims = yLims, N = N, newFigure = 0, plotYLabel = 1, colorBar = 0)
		elif i == len(seqDimensionList) - 1:
			thisPlot = plot1DMultiLine(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, yLims = yLims, N = N, newFigure = 0, plotYLabel = 0, colorBar = colorBar)	
		else:
			thisPlot = plot1DMultiLine(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, yLims = yLims, N = N, newFigure = 0, plotYLabel = 0, colorBar = 0)	
	
	if whatToPlot == 'RR':
		suptitle('Reward Rate')
	elif whatToPlot == 'FC':
		suptitle('Fraction Correct')
		ylim([.5,1])
	elif whatToPlot == 'RT':
		suptitle('Reaction Time')
	
	if colorBar == 1:
		subplots_adjust(bottom=0.12, right=0.85, top=0.8,left=.065)
	else:
		subplots_adjust(bottom=0.12, top=0.8)
	
	if saveFig != 0:
		savefig('/Users/Nick/Desktop/fig1.eps')
	
	return

################################################################################
# This function plots creates a multiline plot:
def plot1DMultiLine(sliceDict, whatToPlot, saveResultDir = './savedResults', whichRun = 0, tDel = 2000, tPen = 0, tND = 350, newFigure = 1, quickName = -1, N = 5, colorBar = 1, titleString = -1, yLims = -1, plotYLabel = 1, color = [], saveFig=0):
	from numpy import array, linspace, inf
	from pylab import figure, subplots_adjust, cm, flipud, pcolor, colorbar, hold, savefig, ylim, close
	import copy
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')
	if whatToPlot == 'All':
		plot1DMultiLine(copy.copy(sliceDict), 'RR', saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, newFigure = newFigure, quickName = quickName, N = N, colorBar = colorBar, titleString = titleString, yLims = yLims, plotYLabel = plotYLabel, color = color)
		plot1DMultiLine(copy.copy(sliceDict), 'RT', saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, newFigure = newFigure, quickName = quickName, N = N, colorBar = colorBar, titleString = titleString, yLims = yLims, plotYLabel = plotYLabel, color = color)
		plot1DMultiLine(copy.copy(sliceDict), 'FC', saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, newFigure = newFigure, quickName = quickName, N = N, colorBar = colorBar, titleString = titleString, yLims = yLims, plotYLabel = plotYLabel, color = color)
		return

	seqDimensionTuple = sliceDict['MultiLineVar']
	if isinstance(seqDimensionTuple, str):
		seqDimension = seqDimensionTuple
		settings, FD, numberOfJobs, gitVersion =  getSettings(quickName, saveResultDir, whichRun=whichRun)
		vals = settings[seqDimension]
		seqDimensionList = linspace(min(vals), max(vals), N)
	else:
		seqDimension = seqDimensionTuple[0]
		seqDimensionList = seqDimensionTuple[1]
		
	if color == []:
		seqDimensionListArray = array(seqDimensionList, dtype=float)
		colorMatrix=cm.autumn_r((seqDimensionListArray-min(seqDimensionListArray))/(max(seqDimensionListArray)-min(seqDimensionListArray)))
	if colorBar: pcolor(array([[min(seqDimensionList),max(seqDimensionList)]]),cmap=cm.autumn_r,visible=False)

	del sliceDict['MultiLineVar']
	if newFigure:
		figure()
	minY = inf
	maxY = -inf
	figure(99)
	for i in range(len(seqDimensionList)):
		sliceDict[seqDimension] = seqDimensionList[i]
		thisPlot = plot1D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = -1, newFigure = 0, plotYLabel = plotYLabel, yLims = -1, color = colorMatrix[i])
		if min(thisPlot[0].get_ydata()) < minY:
			minY = min(thisPlot[0].get_ydata())
		if max(thisPlot[0].get_ydata()) > maxY:
			maxY = max(thisPlot[0].get_ydata())
			
	if yLims == -1:
		yLims = (minY, maxY)
	
	figure(1)		
	for i in range(len(seqDimensionList)):
		sliceDict[seqDimension] = seqDimensionList[i]
		if titleString == -1:
			titleString = ''
		thisPlot = plot1D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, yLims = yLims, newFigure = 0, plotYLabel = plotYLabel, color = colorMatrix[i])	
		
	if colorBar: 
		cb = colorbar()
		cb.set_label('Color variable: ' + seqDimension)
		
	if whatToPlot == 'FC':
		ylim([.5,1])
		
	if saveFig != 0:
		savefig('/Users/Nick/Desktop/fig1.eps')
	close(99)
	
	
	return yLims

################################################################################
# This function plots a sequence  of 1-D plots:
def plot1DSeq(sliceDict, whatToPlot, saveResultDir = './savedResults', whichRun = 0, tDel = 2000, tPen = 0, tND = 350, newFigure = 1, quickName = -1, seqLength = 4):
	from numpy import array, linspace, inf
	from pylab import figure, subplot, suptitle, subplots_adjust, ylim
	import copy
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	seqDimensionTuple = sliceDict['SeqVar']
	if isinstance(seqDimensionTuple, str):
		seqDimension = seqDimensionTuple
		settings, FD, numberOfJobs, gitVersion =  getSettings(quickName, saveResultDir, whichRun=whichRun)
		vals = settings[seqDimension]
		seqDimensionList = linspace(min(vals), max(vals), seqLength)
	else:
		seqDimension = seqDimensionTuple[0]
		seqDimensionList = seqDimensionTuple[1]

	del sliceDict['SeqVar']
	if newFigure:
		figure(num=None,figsize=(4*seqLength, 4))
	minY = inf
	maxY = -inf
	for i in range(len(seqDimensionList)):
		subplot(1,len(seqDimensionList),i+1)
		sliceDict[seqDimension] = seqDimensionList[i]
		titleString = seqDimension + '=' + '%-5.3f' % seqDimensionList[i]
		if i==0:
			thisPlot = plot1D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, newFigure = 0, plotYLabel = 1, yLims = -1)
		else:
			thisPlot = plot1D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, newFigure = 0, plotYLabel = 0, yLims = -1)	
		if min(thisPlot[0].get_ydata()) < minY:
			minY = min(thisPlot[0].get_ydata())
		if max(thisPlot[0].get_ydata()) > maxY:
			maxY = max(thisPlot[0].get_ydata())
			
	yLims = (minY, maxY)
			
	for i in range(len(seqDimensionList)):
		subplot(1,len(seqDimensionList),i+1)
		sliceDict[seqDimension] = seqDimensionList[i]
		titleString = seqDimension + '=' + '%-5.3f' % seqDimensionList[i]
		if i==0:
			thisPlot = plot1D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, yLims = yLims, newFigure = 0, plotYLabel = 1)
		else:
			thisPlot = plot1D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, quickName = quickName, titleString = titleString, yLims = yLims, newFigure = 0, plotYLabel = 0)	
		if whatToPlot == 'FC':
			ylim([.5,1])
	
	if whatToPlot == 'RR':
		suptitle('Reward Rate')
	elif whatToPlot == 'FC':
		suptitle('Fraction Correct')
	elif whatToPlot == 'RT':
		suptitle('Reaction Time')
		
	subplots_adjust(bottom=0.12, right=0.97, top=0.8,left=.065)
		
	return

################################################################################
# This function plots a sequence  of 2-D slices:
def plot2DSeq(sliceDict, whatToPlot, saveResultDir = './savedResults', whichRun = 0, tDel = 2000, tPen = 0, tND = 350, newFigure = 1, colorArray = [], N = 20, quickName = -1, seqLength = 4,  colorBar = 1):
	from numpy import array, linspace, inf
	from pylab import figure, subplot, colorbar, suptitle, subplots_adjust
	import copy
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	seqDimensionTuple = sliceDict['SeqVar']
	if isinstance(seqDimensionTuple, str):
		seqDimension = seqDimensionTuple
		settings, FD, numberOfJobs, gitVersion =  getSettings(quickName, saveResultDir, whichRun=whichRun)
		vals = settings[seqDimension]
		seqDimensionList = linspace(min(vals), max(vals), seqLength)
	else:
		seqDimension = seqDimensionTuple[1]
		seqDimensionList = seqDimensionTuple[0]

	del sliceDict['SeqVar']
	if newFigure:
		figure(num=None,figsize=(4*seqLength, 4))
	minZ = inf
	maxZ = -inf
	for i in range(len(seqDimensionList)):
		subplot(1,len(seqDimensionList),i+1)
		sliceDict[seqDimension] = seqDimensionList[i]
		titleString = seqDimension + '=' + '%-5.3f' % seqDimensionList[i]
		if i+1 == len(seqDimensionList):
			thisPlot = plot2D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, N = N, quickName = quickName, titleString = titleString, colorBar = 0, newFigure = 0, plotYLabel = 0)
		elif i==0:
			thisPlot = plot2D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, N = N, quickName = quickName, titleString = titleString, colorBar = 0, newFigure = 0, plotYLabel = 0)
		else:
			thisPlot = plot2D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, N = N, quickName = quickName, titleString = titleString, colorBar = 0, newFigure = 0, plotYLabel = 0)	
		if min(thisPlot.levels) < minZ:
			minZ = min(thisPlot.levels)
		if max(thisPlot.levels) > maxZ:
			maxZ = max(thisPlot.levels)
	
	colorArray = linspace(minZ, maxZ, N)
			
	for i in range(len(seqDimensionList)):
		subplot(1,len(seqDimensionList),i+1)
		sliceDict[seqDimension] = seqDimensionList[i]
		titleString = seqDimension + '=' + '%-5.3f' % seqDimensionList[i]
		if i+1 == len(seqDimensionList):
			thisPlot = plot2D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, colorArray = colorArray, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, N = N, quickName = quickName, titleString = titleString, colorBar = 1, newFigure = 0, plotYLabel = 0)
		elif i==0:
			thisPlot = plot2D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, colorArray = colorArray, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, N = N, quickName = quickName, titleString = titleString, colorBar = 0, newFigure = 0, plotYLabel = 1)
		else:
			thisPlot = plot2D(copy.copy(sliceDict), whatToPlot,saveResultDir = saveResultDir, colorArray = colorArray, whichRun = whichRun, tDel = tDel, tPen = tPen, tND = tND, N = N, quickName = quickName, titleString = titleString, colorBar = 0, newFigure = 0, plotYLabel = 0)	
	
	if whatToPlot == 'RR':
		suptitle('Reward Rate')
	elif whatToPlot == 'FC':
		suptitle('Fraction Correct')	
	elif whatToPlot == 'RT':
		suptitle('Reaction Time')
		
	subplots_adjust(bottom=0.12, right=0.85, top=0.8,left=.065)
		
	return

################################################################################
# This function plots a 2-D slice:
def plot2D( sliceDict, whatToPlot,saveResultDir = './savedResults', whichRun = 0, tDel = 2000, tPen = 0, tND = 350, newFigure = 1, colorArray = [], N = 20, quickName = -1, colorBar = 1, plotYLabel = 1, titleString = -1):
	from numpy import transpose, shape, squeeze, array
	import pylab as pl
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	# Get data:
	crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion =  getDataAndSettings(quickName, saveResultDir, whichRun)
	crossTimeData += tND 
	
	# Record variable to plot, and then strip input dictionary of that variable:
	xDimension = sliceDict['XVar']
	yDimension = sliceDict['YVar']
	del sliceDict['XVar']
	del sliceDict['YVar']
	
	# Reorder dimension list and cube to put plotting variable first:
	permuteList = range(len(dims))
	whereIsXDim = dims.index(xDimension)
	whereIsYDim = dims.index(yDimension)
	dims[1], dims[whereIsXDim] = dims[whereIsXDim], dims[1]
	dims[0], dims[whereIsYDim] = dims[whereIsYDim], dims[0]
	permuteList[1], permuteList[whereIsXDim] = permuteList[whereIsXDim], permuteList[1]
	permuteList[0], permuteList[whereIsYDim] = permuteList[whereIsYDim], permuteList[0]
	crossTimeData = transpose(crossTimeData,permuteList)
	resultData = transpose(resultData,permuteList)
	
	# Collapse all non-constant dimensions:
	crossDims = dims[:]
	resultDims = dims[:]
	for collapseDim in iter(sliceDict):
		crossTimeData, resultData, crossDims = reduce1D(crossTimeData, resultData, crossDims, collapseDim, settings[collapseDim], sliceDict[collapseDim], tDel = tDel, tPen = tPen, tND = tND)
	crossTimeSlice = squeeze(crossTimeData)
	resultSlice = squeeze(resultData)
	
	# Create x-axis values, and plot:
	xVals = settings[xDimension]
	yVals = settings[yDimension]
	if whatToPlot == 'RR':
		depVar = 1000*resultSlice/(crossTimeSlice + tND + tDel + (1-resultSlice)*tPen)
		heightLabel = 'Reward Rate'
	elif whatToPlot == 'RT':
		depVar = crossTimeSlice
		heightLabel = 'Reaction Time'
	elif whatToPlot == 'FC':
		depVar = resultSlice
		heightLabel = 'Fraction Correct'
	else: print ' Unrecognized plot option ' + whatToPlot
	
	if newFigure:
		pl.figure()
	if colorArray == []:
		myPlot = pl.contourf(xVals,yVals,depVar,N)
	else:
		myPlot = pl.contourf(xVals,yVals,depVar,N,levels = colorArray)
	pl.xlabel(xDimension)
	if plotYLabel:
		pl.ylabel(yDimension)
	if titleString == -1:
		pl.title(heightLabel)
	else:
		pl.title(titleString)
	if colorBar:
		pl.colorbar()
	return myPlot

################################################################################
# This function plots a 1-D slice:
def plot1D( sliceDict, whatToPlot,saveResultDir = './savedResults', whichRun = 0, tDel = 2000, tPen = 0, tND = 350, quickName = -1, titleString = -1, newFigure = 1, plotYLabel = 1, yLims = -1, color = -1, saveFig=0,):
	from numpy import transpose, shape, squeeze, ndarray, array, mean, exp
	import pylab as pl
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	# Get data:
	crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion =  getDataAndSettings(quickName, saveResultDir, whichRun)
	crossTimeData += tND 
	
	# Record variable to plot, and then strip input dictionary of that variable:
	xDimension = sliceDict['XVar']
	del sliceDict['XVar']
	
	# Reorder dimension list and cube to put plotting variable first:
	permuteList = range(len(dims))
	whereIsXDim = dims.index(xDimension)
	dims[0], dims[whereIsXDim] = dims[whereIsXDim], dims[0]
	permuteList[0], permuteList[whereIsXDim] = permuteList[whereIsXDim], permuteList[0]
	crossTimeData = transpose(crossTimeData,permuteList)
	resultData = transpose(resultData,permuteList)
    
	# Reorder dimension list and cube to put theta variable last:
	if 'theta' in sliceDict:
		permuteList = range(len(dims))
		whereIsXDim = dims.index('theta')
		dims[-1], dims[whereIsXDim] = dims[whereIsXDim], dims[-1]
		permuteList[-1], permuteList[whereIsXDim] = permuteList[whereIsXDim], permuteList[-1]
		crossTimeData = transpose(crossTimeData,permuteList)
		resultData = transpose(resultData,permuteList)
	
	# Collapse all non-constant dimensions:
	crossDims = dims[:]
	resultDims = dims[:]
	for collapseDim in iter(sliceDict):
		crossTimeData, resultData, crossDims = reduce1D(crossTimeData, resultData, crossDims, collapseDim, settings[collapseDim], sliceDict[collapseDim], tDel = tDel, tPen = tPen, tND = tND)
	crossTimeSlice = squeeze(crossTimeData)
	resultSlice = squeeze(resultData)
	
	# Create x-axis values, and plot:
	xVals = settings[xDimension]
	if whatToPlot == 'RR':
		depVar = 1000*resultSlice/(crossTimeSlice + tND + tDel + (1-resultSlice)*tPen)
		yAxisLabel = 'Reward Rate'
	elif whatToPlot == 'RT':
		depVar = crossTimeSlice
		yAxisLabel = 'Reaction Time'
	elif whatToPlot == 'FC':
		depVar = resultSlice
		yAxisLabel = 'Fraction Correct'
		yLims = (0,1)
	else: print ' Unrecognized plot option ' + whatToPlot
	
	if newFigure:
		pl.figure()
	if not(isinstance(color,ndarray)):
		myPlot = pl.plot(xVals,depVar)
	else:
		myPlot = pl.plot(xVals,depVar, color = color)
	pl.xlim((min(xVals),max(xVals)))
	
	# Add RT and ER physiological data is appropriate:
	if FD == 1:
		inputSet = 'FD'
	else:
		inputSet = 'FR'
	if (whatToPlot == 'FC') and ((xDimension == 'xMean') or (xDimension == 'C') or (xDimension == 'COn') or (xDimension == 'CPre')):
		CData=array([0,3.2,6.4,12.8,25.6,51.2])
		FC, RT = getRoitmanPsyChr(inputSet)
		pl.plot(CData,FC,'o')
	elif (whatToPlot == 'RT') and ((xDimension == 'xMean') or (xDimension == 'C') or (xDimension == 'COn') or (xDimension == 'CPre')) and (FD == 0):
		CData=array([0,3.2,6.4,12.8,25.6,51.2])
		FC, RT = getRoitmanPsyChr(inputSet)
		pl.plot(CData, RT,'o')
	
	if yLims != -1:
		pl.ylim(yLims[0], yLims[1])
	pl.xlabel(xDimension)
	if plotYLabel:
		pl.ylabel(yAxisLabel)
	if titleString != -1:
		pl.title(titleString)
	
	if saveFig != 0:
		pl.savefig('/Users/Nick/Desktop/fig1.eps')
	
	return myPlot#xVals,depVar
	
################################################################################
# This function reduces the dimension of a cube by 1, along a given slice:
def reduce1D(crossTimeCube, resultCube, dims, varToReduce, vals, sliceVal, tDel = 2000, tPen = 0, tND = 350):

	if varToReduce == 'theta' and sliceVal == 'Optimize':
		reduceTuple = reduceThetaOptimize(crossTimeCube, resultCube, dims, tDel = tDel, tPen = tPen, tND = tND)
		return reduceTuple
	if isinstance(sliceVal,list) and sliceVal[0] == 'Marginalize':
		reduceTuple = reduceMarginalize(crossTimeCube, resultCube, dims, varToReduce, vals, sliceVal[1],sliceVal[2])
		return reduceTuple
	else:
		for i in range(len(vals)-1):
			if vals[i] <= sliceVal and sliceVal < vals[i+1]:
				lInd=i
				rInd=i+1
				break
		indexListL = [slice(None,None)]*len(dims)
		indexListR = indexListL[:]
		indToSet = dims.index(varToReduce)
		if vals[-1] != sliceVal:
			indexListL[indToSet] = lInd
			indexListR[indToSet] = rInd
			crossTimeCubeL = crossTimeCube[tuple(indexListL)]
			crossTimeCubeR = crossTimeCube[tuple(indexListR)]
			crossTimeCubeReduce = (crossTimeCubeL*float(vals[rInd] - sliceVal) + crossTimeCubeR*float(sliceVal - vals[lInd]))/float(vals[rInd] - vals[lInd])
			resultCubeL = resultCube[tuple(indexListL)]
			resultCubeR = resultCube[tuple(indexListR)]
			resultCubeReduce = (resultCubeL*float(vals[rInd] - sliceVal) + resultCubeR*float(sliceVal - vals[lInd]))/float(vals[rInd] - vals[lInd])
		else:
			indexListR[indToSet] = -1
			crossTimeCubeReduce = crossTimeCube[tuple(indexListR)]
			resultCubeReduce = resultCube[tuple(indexListR)]
		dims.remove(varToReduce)
		return (crossTimeCubeReduce, resultCubeReduce, dims)

################################################################################
# This function reduces a dimension by marginalizing:
def reduceMarginalize(crossTimeCube, resultCube, dims, varToReduce, vals, distribution, distributionSettings):

	# Import necessary packages:
	import copy
	import numpy as np

	# Pick Distribution:
	numOfValues = len(vals)
	vals = np.array(vals)
	if distribution == 'Normal':
		delta = vals[1]-vals[0]
		vals = vals - delta*1./2
		L = np.append(vals[1]-10000, vals[1:])
		R = np.append(vals[1:], vals[-1]+10000)
		probDists = [0]*len(L)
		for j in range(len(L)):
			probDists[j] = intErfAB(L[j], R[j], mu=distributionSettings[0], sigma=distributionSettings[1])
	elif distribution == 'Delta':
		probDists = np.zeros(numOfValues)
		deltaIndex = np.argmin(np.abs(vals-distributionSettings))
		probDists[deltaIndex] = 1
	elif distribution == 'Uniform':
		probDists = np.zeros(numOfValues)        
		deltaIndexL = np.argmin(np.abs(vals-distributionSettings[0]))
		deltaIndexR = np.argmin(np.abs(vals-distributionSettings[1]))
		numNonZero = deltaIndexR - deltaIndexL + 1
		probDists[deltaIndexL:(deltaIndexR + 1)] = 1./numNonZero
	else:
		print 'Not Supported: ' + distribution
		import sys
		sys.exit()
        
    # Double-check the distribution:
	if np.sum(probDists) < .999:
		print 'pmf doent sum to 1'
		print np.sum(probDists)
		import sys
		sys.exit()
		
	# Marginalize across variable:
	marginalVal = vals[0]
	marginalCrossTimeCube, marginalResultCube, dimsMarginal = reduce1D(crossTimeCube, resultCube, copy.copy(dims), varToReduce, vals, marginalVal)
	marginalCrossTimeCube *= probDists[0]
	marginalResultCube *= probDists[0]
	for j in range(1,numOfValues):
		marginalVal = vals[j]
		crossTimeCubeTemp, resultCubeTemp, dimsMarginalTemp = reduce1D(crossTimeCube, resultCube, copy.copy(dims), varToReduce, vals, marginalVal)
		marginalCrossTimeCube += probDists[j]*crossTimeCubeTemp
		marginalResultCube += probDists[j]*resultCubeTemp

	return (marginalCrossTimeCube, marginalResultCube, dimsMarginal)
	
################################################################################
# This function reduces the Theta dimension through RR optimization:
def reduceThetaOptimize(crossTimeCube, resultCube, dims, tDel = 2000, tPen = 0, tND = 350):
	import numpy as np
	import itertools

	# Create RR Cube:
	RRCube = resultCube/(crossTimeCube + tND + tDel + (1-resultCube)*tPen)
	
	# Cut down input Arrays by the theta dimension:
	indexList = [slice(None,None)]*len(dims)
	indToSet = dims.index('theta')
	indexList[indToSet] = 0
	crossTimeCubeReduce = crossTimeCube[tuple(indexList)]
	resultCubeReduce = resultCube[tuple(indexList)]
	
	# Pick Optimal entries:
	indexList = [slice(None,None)]*len(dims)
	indToSet = dims.index('theta')
	iterList = [range(dimLength) for dimLength in np.shape(RRCube)]
	iterList[indToSet] = [0]
	for indexTuple in itertools.product(*iterList):
		indexList = list(indexTuple)
		indexListOut = indexList[:]
		indexList[indToSet] = slice(None,None)
		tmp = indexListOut.pop(indToSet)
		currThetaList = RRCube[tuple(indexList)]
		indexList[indToSet] = currThetaList.argmax()
		crossTimeCubeReduce[tuple(indexListOut)] = crossTimeCube[tuple(indexList)]
		resultCubeReduce[tuple(indexListOut)] = resultCube[tuple(indexList)]
		
	dims.remove('theta')
	return (crossTimeCubeReduce, resultCubeReduce, dims)
	

################################################################################
# This function lists the job names that are available:
def listNames(saveResultDir = './savedResults', N=10):

	import operator
	nameDict = quickNameIDDictionary(saveResultDir,includeRepeats = 0)
	nameTimeList = []
	for item in nameDict:
		nameTimeList.append((item, nameDict[item][0][1]))
	nameTimeListSorted = sorted(nameTimeList, key=operator.itemgetter(1),reverse=True)
	print ' Available job names:'
	counter = 0
	for name in nameTimeListSorted:
		counter += 1
		if counter <= N:
			print '   ' + name[0] + ', (' + str(getTrials(quickName=name[0])) + ')'
	return

################################################################################
# This function prints out a nicely formatted settings string:
def printSettings(quickName = -1, saveResultDir = './savedResults', whichRun = 0):

	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	printString = getSettingsString(quickName, saveResultDir = saveResultDir, whichRun = whichRun)
	print printString
	return

################################################################################
# This function gets the settings string from a file:
def getSettingsString(quickName = -1, saveResultDir = './savedResults', whichRun = 0):
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	settings, FD, numberOfJobs, gitVersion = getSettings(quickName, saveResultDir, whichRun=whichRun)
	params = settings.keys()
	constParams = []
	varParams = []
	for parameter in params:
		if len(settings[parameter])>1: varParams.append(parameter)
		else: constParams.append(parameter)
	constParams.sort
	varParams.sort
	settingsString = ' Job "quickName": ' + quickName + '\n'
	if FD:
		settingsString += ' Fixed-Duration protocol (FD)\n'
	else:
		settingsString += ' Reaction-Time protocol (RT)\n'
	settingsString += ' Parameter Settings:\n'
	totalLength = 1
	for parameter in constParams:
		thisSetting = settings[parameter]
		settingsString += '   %6s: %10.2f\n' % (parameter, min(thisSetting))
		totalLength *= len(thisSetting)
	for parameter in varParams:
		thisSetting = settings[parameter]
		settingsString += '   %6s: %10.2f %5.2f %3d\n' % (parameter,min(thisSetting),max(thisSetting),len(thisSetting))
		totalLength *= len(thisSetting)

	settingsString += ' Drift-Diffusion Software version:  %-5s\n' % gitVersion		
	settingsString += ' Number of Parameter Space Points: %-5d\n' % totalLength
	settingsString += ' Number of Simulations per Point:  %-5d\n' % numberOfJobs[1]
	settingsString += ' Total number of Simulations:      %-5d' % (totalLength*numberOfJobs[1])
	
	return settingsString


################################################################################
# This function gets the name of a file, given an ID:
def getFileString(ID, typeOfFile,  saveResultDir = './savedResults'):

	resultDict = IDquickNameDictionary(saveResultDir)
	quickName = resultDict[ID]
	fileName = quickName + '_' + ID + '.' + typeOfFile
	return fileName

################################################################################
# This function grabs the data for a given quickName:
def getData(quickName = -1, saveResultDir = './savedResults', whichRun = 0):
	import pickle
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')
	
	ID = quickNameToID(quickName=quickName, saveResultDir=saveResultDir, whichRun=whichRun)
	fileName = getFileString(ID,'dat', saveResultDir)
	fIn = open(saveResultDir + '/' + fileName,'r')
	resultTuple = pickle.load(fIn)
	return resultTuple

################################################################################
# This function grabs the settings for a given quickName:
def getSettings(quickName = -1, saveResultDir = './savedResults', whichRun = 0):
	import pickle
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')
		
	ID = quickNameToID(quickName=quickName, saveResultDir=saveResultDir, whichRun = whichRun)
	fileName = getFileString(ID,'settings', saveResultDir)
	fIn = open(saveResultDir + '/' + fileName,'r')
	resultTuple = pickle.load(fIn)
	return resultTuple

################################################################################
# This function grabs the results and settings for a given quickName:
def getDataAndSettings(quickName = -1, saveResultDir = './savedResults', whichRun = 0):
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')
	
	crossTimeData, resultData, dims = getData(quickName=quickName, saveResultDir=saveResultDir, whichRun=whichRun)
	settings, FD, numberOfJobs, gitVersion = getSettings(quickName, saveResultDir, whichRun=whichRun)
	return (crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion)
	
################################################################################
# This function creates a dictionary between file "ID's" and quickNames
def IDquickNameDictionary(saveResultDir = './savedResults'):
	import pickle, os, operator
	
	resultDict = {}
	for root, dirs, files in os.walk(saveResultDir):
		for name in files:
			if len(name.split('.')) > 1 and len(name.split('.')) < 3:
				quickNameAndID,suffix = name.split('.')
				if suffix == 'settings':
					quickName, ID = quickNameAndID.split('_')
					resultDict[ID] = quickName
	return resultDict

################################################################################
# This function creates a dictionary between "quickNames" and file ID's
def quickNameIDDictionary(saveResultDir = './savedResults',includeRepeats = 0):
	import pickle, os, operator
	
	resultDict = {}
	for root, dirs, files in os.walk(os.path.abspath(saveResultDir)):
		for name in files:
			if len(name.split('.')) > 1 and len(name.split('.')) < 3:
				quickNameAndID,suffix = name.split('.')
				if suffix == 'settings':
					st = os.stat(os.path.join(root, name))
					quickName, ID = quickNameAndID.split('_')
					IDTime = st[8]
					fIn = open(os.path.join(root, name),'r')
					settingTuple = pickle.load(fIn)

					if not(resultDict.has_key(quickName)):
						resultDict[quickName] = [(ID, IDTime)]
					else:
						tempList = resultDict[quickName]
						tempList.append((ID, IDTime))
						tempListSorted = sorted(tempList, key=operator.itemgetter(1))
						if includeRepeats:
							resultDict[quickName] = tempListSorted
						else:
							resultDict[quickName] = [tempListSorted[-1]]
	return resultDict

################################################################################
# This function grabs the ID for a given quickName:
def quickNameToID(quickName = -1, saveResultDir = './savedResults', whichRun = 0):
	import operator
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	currentDict = quickNameIDDictionary(saveResultDir=saveResultDir, includeRepeats = 1)
	try: listOfIDTimeTuple = currentDict[quickName]
	except KeyError: 
		print '  Job "' + quickName + '" not found.'
		print '  Available jobs:'
		for i in currentDict.keys(): print '    ' + i
		raise
	listOfID = map(operator.itemgetter(0), listOfIDTimeTuple)
	return listOfID[whichRun-1]
		
################################################################################
# This function gets the most recent quickname:
def getLastQuickName(saveResultDir = './savedResults'):
	import operator
	
	d = quickNameIDDictionary(saveResultDir=saveResultDir)
	d2 = IDquickNameDictionary(saveResultDir=saveResultDir) 
	myIndex = [d[key][0] for key in iter(d)]
	myIndexSorted = sorted(myIndex, key=operator.itemgetter(1))
	IDName = myIndexSorted[-1][0]
	lastQuickName = d2[IDName]
	return lastQuickName
	
################################################################################
# This function gets the number of trials for a given quickname:
def getTrials(quickName = -1, saveResultDir = './savedResults'):
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')
		
	d = quickNameIDDictionary(saveResultDir = saveResultDir, includeRepeats = 1)
	return len(d[quickName])
	
################################################################################
# This function loads the Roitman Data Set, psychr:
def getRoitmanPsyChr(inputSet):

	import os
	dataSetDir = '/Users/Nick/Desktop/currentProjects/DDMCubeTeragrid/'

	if inputSet == 'FR':
		psyChrFileName = 'Roitman_data_psychoCronoData.dat'
	elif inputSet == 'FD':
		psyChrFileName = 'Roitman_data_psychoCronoDataFD.dat'
	
	f = open(os.path.join(dataSetDir,psyChrFileName), 'r')
	FC = [0]*6
	RT = [0]*6
	counter = -1
	for line in f:
		if line[0] != 'c' and line != '\n':
			counter = counter + 1
			coh,FC[counter],RT[counter] = line.split('	')
	f.close()
	
	for i in range(len(FC)):
		FC[i] = float(FC[i])
		RT[i] = float(RT[i])
	
	return FC,RT
	
################################################################################
# This function loads the Roitman Data Set, RTCurve:
def getRoitmanRTCurve():
	RTCurveFileName = 'Roitman_data_RTCurve_Coh_6.4.dat'
	
	return
	
################################################################################
# This function exports a 1-D slice:
def export1D( sliceDict, whatToPlot,saveResultDir = './savedResults', whichRun = 0, tDel = 2000, tPen = 0, tND = 350, quickName = -1, titleString = -1, plotYLabel = 1):
	from numpy import transpose, shape, squeeze, ndarray, array, mean, exp
	import pylab as pl, copy
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	# Get data:
	crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion =  getDataAndSettings(quickName=quickName, saveResultDir=saveResultDir, whichRun=whichRun)
	crossTimeData += tND 
	
	# Record variable to plot, and then strip input dictionary of that variable:
	xDimension = sliceDict['XVar']
	del sliceDict['XVar']
	
	# Reorder dimension list and cube to put plotting variable first:
	permuteList = range(len(dims))
	whereIsXDim = dims.index(xDimension)
	dims[0], dims[whereIsXDim] = dims[whereIsXDim], dims[0]
	permuteList[0], permuteList[whereIsXDim] = permuteList[whereIsXDim], permuteList[0]
	crossTimeData = transpose(crossTimeData,permuteList)
	resultData = transpose(resultData,permuteList)
	
	# Collapse all non-constant dimensions:
	crossDims = dims[:]
	resultDims = dims[:]
	for collapseDim in iter(sliceDict):
		crossTimeData, resultData, crossDims = reduce1D(crossTimeData, resultData, crossDims, collapseDim, settings[collapseDim], sliceDict[collapseDim], tDel = tDel, tPen = tPen, tND = tND)
	crossTimeSlice = squeeze(crossTimeData)
	resultSlice = squeeze(resultData)
	
	# Create x-axis values, and plot data:
	xVals = settings[xDimension]
	if whatToPlot == 'RR':
		depVar = 1000*resultSlice/(crossTimeSlice + 0 + tDel + (1-resultSlice)*tPen)
	elif whatToPlot == 'RT':
		depVar = crossTimeSlice
	elif whatToPlot == 'FC':
		depVar = resultSlice
	else: print ' Unrecognized plot option ' + whatToPlot

	return xVals, depVar

################################################################################
# This function exports a 1-D slice:
def export1Elem( sliceDict, whatToReturn = 'both',saveResultDir = './savedResults', whichRun = 0, quickName = -1,tND = 350):
	from numpy import transpose, shape, squeeze, ndarray, array, mean, exp
	import pylab as pl, copy
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	# Get data:
	crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion =  getDataAndSettings(quickName, saveResultDir, whichRun)
	
	# Reorder dimension list and cube:
	permuteList = range(len(dims))
	crossTimeData = transpose(crossTimeData,permuteList)
	resultData = transpose(resultData,permuteList)
	
	# Collapse all non-constant dimensions:
	crossDims = dims[:]
	resultDims = dims[:]
	for collapseDim in iter(sliceDict):
		crossTimeData, resultData, crossDims = reduce1D(crossTimeData, resultData, crossDims, collapseDim, settings[collapseDim], sliceDict[collapseDim])
	crossTimeSlice = squeeze(crossTimeData)
	resultSlice = squeeze(resultData)


	if whatToReturn == 'both':
		return (crossTimeSlice.tolist() + tND, resultSlice.tolist())
	elif whatToReturn == 'FC':
		return resultSlice.tolist()
	elif whatToReturn == 'RT':
		return crossTimeSlice.tolist() + tND
	else:
		return -1

################################################################################
# This function plots a histogram:
def histPlot( sliceDict,saveResultDir = './savedResults', whichRun = 0, quickName = -1,tND = 350, bins=20, normed=True, plotsOn = 1,center=1, CorI = 'Both',newFigure=True):
	from numpy import transpose, shape, squeeze, ndarray, array, mean, exp, atleast_2d, nonzero, mean

	import pylab as pl
	import pickle
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	# Get hash value:
	hashInd = int(export1Elem( sliceDict, whatToReturn = 'FC',saveResultDir = './savedResults', whichRun = whichRun, quickName = quickName,tND = tND))
	
	# Load actual data:
	ID = quickNameToID(quickName, saveResultDir, whichRun=whichRun)
	fileName = getFileString(ID,'dat', saveResultDir)
	first,second = fileName.split('.')
	fileName = first + 'Hash' + '.' +second
	fIn = open(saveResultDir + '/' + fileName,'r')
	resultTuple = pickle.load(fIn)
	myRTData = transpose(atleast_2d(resultTuple[0][hashInd]))
	if CorI == 'Both':
		myRTDataPlot = myRTData + tND
	elif (CorI == 'C') or (CorI == 'I'):
		myFCData = transpose(atleast_2d(resultTuple[1][hashInd]))
		if CorI == 'C':
			targetValue = 1
		else:
			targetValue = 0
		targetIndices = nonzero(myFCData==targetValue)
		myRTDataPlot = myRTData[targetIndices]  + tND
	else:
		print 'Unrecognized option for CorI: ' + CorI
		from sys import exit
		exit(1)

	# Center if desired:
	if center == 1:
		myRTDataPlot = myRTDataPlot - mean(myRTDataPlot)
	
	# Plot the histogram:
	if plotsOn:
		if newFigure: pl.figure()
		else: pl.figure(1).clear()
		histOut = pl.hist(myRTDataPlot,bins=bins,normed=normed)
	else:
		pl.figure(99)
		histOut = pl.hist(myRTDataPlot,bins=bins,normed=normed)
		pl.close(99)
	
	return histOut, myRTDataPlot
	
################################################################################
# This function plots a two histograms, for correct and incorrect results:
def histPlotCICompare(sliceDict,saveResultDir = './savedResults', whichRun = 0, quickName = -1,tND = 350, bins=20, normed=True, plotsOn = 1,center=1,newFigure=True):

	# Import necessary functions:
	import pylab as pl
	from numpy import mean

	# Gather data:
	histOutC, myCorrData = histPlot(sliceDict,saveResultDir=saveResultDir,whichRun=whichRun,quickName=quickName,tND=tND,bins=bins,normed=normed,plotsOn=0,center=center,CorI='C')
	histOutI, myICorrData = histPlot(sliceDict,saveResultDir=saveResultDir,whichRun=whichRun,quickName=quickName,tND=tND,bins=bins,normed=normed,plotsOn=0,center=center,CorI='I')
	
	# Create plotting data, correct:
	corrX = histOutC[1][0:-1]
	widthC = corrX[1]-corrX[0]
	corrY = histOutC[0]
	
	# Create plotting data, incorrect:
	iCorrX = histOutI[1][0:-1]
	widthI = iCorrX[1]-iCorrX[0]
	iCorrY = -1*histOutI[0]
	
	# Make the double-plot:
	if newFigure: pl.figure()
	pl.bar(corrX,corrY,width=widthC,color='green')
	pl.bar(iCorrX,iCorrY,width=widthI,color='red')
	
	# Print mean of each dist:
	print 'Mean, correct: ' + str(mean(myCorrData))
	print 'Mean, incorrect: ' + str(mean(myICorrData))
	
	return

################################################################################
# This function plots a multi-histogram plot:
def histPlotMultiBar(sliceDict,saveResultDir = './savedResults', whichRun = 0, quickName = -1,tND = 350, bins=20, normed=True, plotsOn = 1,center=1,N=2, newFigure=True):

	from numpy import transpose, shape, squeeze, ndarray, array, mean, exp, atleast_2d, linspace, concatenate
	import pylab as pl
	import pickle
	if quickName == -1:
		quickName = getLastQuickName(saveResultDir = './savedResults')

	# Get data:
	crossTimeData, resultData, dims, settings, FD, numberOfJobs, gitVersion =  getDataAndSettings(quickName, saveResultDir, whichRun)
	
	# Get multiBar stats:
	seqDimensionTuple = sliceDict['MultiBarVar']
	if isinstance(seqDimensionTuple, str):
		seqDimension = seqDimensionTuple
		settings, FD, numberOfJobs, gitVersion =  getSettings(quickName, saveResultDir, whichRun=whichRun)
		vals = settings[seqDimension]
		seqDimensionList = linspace(min(vals), max(vals), N)
	else:
		seqDimension = seqDimensionTuple[0]
		seqDimensionList = seqDimensionTuple[1]
	del sliceDict['MultiBarVar']

	# Get hash value:
	hashInd = [0]*len(seqDimensionList)
	for i in range(len(hashInd)):
		sliceDict[seqDimension] = seqDimensionList[i]
		hashInd[i] = int(export1Elem( sliceDict, whatToReturn = 'FC',saveResultDir = './savedResults', whichRun = whichRun, quickName = quickName,tND = tND))
	
	# Load actual data:
	ID = quickNameToID(quickName, saveResultDir, whichRun=whichRun)
	fileName = getFileString(ID,'dat', saveResultDir)
	first,second = fileName.split('.')
	fileName = first + 'Hash' + '.' +second
	fIn = open(saveResultDir + '/' + fileName,'r')
	resultTuple = pickle.load(fIn)
	
	myRTData = [0]*len(seqDimensionList)
	for i in range(len(hashInd)):
		myRTData[i] = transpose(atleast_2d(resultTuple[0][hashInd[i]]))

		# Center if desired:
		if center == 1:
			myRTData[i] = myRTData[i] - mean(myRTData[i])
	
	# Plot the histogram:
	if newFigure: 
		pl.figure()
	else:
		pl.figure().clear()		
	if plotsOn:

		pl.hist(concatenate(myRTData,axis=1),bins=bins,normed=normed)
	
	return myRTData
	
	
################################################################################
# This function determines error between a slice and data:
def sliceErr(sliceDict, FRorFD='FR', saveResultDir='./savedResults', quickName=-1, whichRun=0, tND=350):
	
	# Import packages
	import copy
	import numpy as np
	import pylab as pl
	
	# Settings:
	FCWeight = .5
	RTWeight = .5
	
	# Set x-axis to coherence:
	sliceDict['XVar'] = 'COn'
	
	# Get  data:
	xValsMonkey = [0,3.2,6.4,12.8,25.6,51.2]
	if FRorFD == 'FR':

		# Simulation:
		xVals, FCData = export1D(copy.copy(sliceDict),'FC',saveResultDir=saveResultDir, quickName=quickName, whichRun=whichRun, tND=tND)
		xVals, RTData = export1D(copy.copy(sliceDict),'RT',saveResultDir=saveResultDir, quickName=quickName, whichRun=whichRun, tND=tND)
		
		# Monkey data:
		FCDataMonkey,RTDataMonkey = getRoitmanPsyChr('FR')
		
		# Interpolate simulated data to hit the correct xValsMonkey points:
		FCDataInterp = np.interp(xValsMonkey, xVals, FCData)
		RTDataInterp = np.interp(xValsMonkey, xVals, RTData)
		
		# Normalize all data sets between 0 and 1, FC:
		FCDataMonkeyNorm = (FCDataMonkey - np.min(FCDataMonkey))/(np.max(FCDataMonkey) - np.min(FCDataMonkey))
		FCDataInterpNorm = (FCDataInterp - np.min(FCDataMonkey))/(np.max(FCDataMonkey) - np.min(FCDataMonkey))
		
		# Normalize all data sets between 0 and 1, RT:
		RTDataMonkeyNorm = (RTDataMonkey - np.min(RTDataMonkey))/(np.max(RTDataMonkey) - np.min(RTDataMonkey))
		RTDataInterpNorm = (RTDataInterp - np.min(RTDataMonkey))/(np.max(RTDataMonkey) - np.min(RTDataMonkey))
		
		# Error calculation:
		errFCAbs = np.sum(np.abs(FCDataInterpNorm - FCDataMonkeyNorm))
		errFC = np.abs(np.sum(FCDataInterpNorm - FCDataMonkeyNorm))
		errRTAbs = np.sum(np.abs(RTDataInterpNorm - RTDataMonkeyNorm))
		errRT = np.abs(np.sum(RTDataInterpNorm - RTDataMonkeyNorm))
		
		# Final error:
		err = FCWeight*(errFCAbs) + RTWeight*(errRTAbs)

	elif FRorFD == 'FD':
		print 'not done yet'
	else:
		print 'unrecoginized option: FRorFD=' + str(FRorFD)
	
	return err

################################################################################
# This function sweeps across noiseSigma and theta, to minimize error:
def minimizeNoiseTheta(theta = 'Free', noiseSigma = 'Free', betaSigma = 0, chopHat = 0, FRorFD='FR', saveResultDir='./savedResults', quickName=-1, whichRun=0, tND=350, saveFig = 0):

	# Import packages:
	import numpy as np
	import copy
	import pylab as pl

	# Settings:
	if betaSigma == 0:
		beta = ['Marginalize','Delta',[0]]
	else:
		beta = ['Marginalize','Normal',[0,betaSigma]]
	uniqueString = 'chopHat: '+str(chopHat)+', beta: '+str(beta)

	# Get settings for theta and noisesigma:
	settingsDict = getSettings(quickName = quickName, saveResultDir = saveResultDir, whichRun = whichRun)[0]
	if theta == 'Free':
		thetaVals = settingsDict['theta']
	else:
		thetaVals = theta
	if noiseSigma == 'Free':
		noiseSigmaVals = settingsDict['noiseSigma']
	else:
		noiseSigmaVals = noiseSigma
	
	# Get error at each point:
	sliceDict = {'chopHat':chopHat,'beta':beta}
	errorMatrix = np.zeros((len(thetaVals),len(noiseSigmaVals)))
	for i in range(len(thetaVals)):
		for j in range(len(noiseSigmaVals)):
			sliceDict['theta'] = thetaVals[i]
			sliceDict['noiseSigma'] = noiseSigmaVals[j]
			errorMatrix[i][j]=sliceErr(copy.copy(sliceDict), FRorFD=FRorFD, saveResultDir=saveResultDir, quickName=quickName, whichRun=whichRun, tND=tND)

	# Generate contour plot
	pl.figure(1).clear()
	pl.figure(1)
	pl.contourf(noiseSigmaVals, thetaVals, errorMatrix, 50)
	pl.title(uniqueString)
	
	# Find minimum element:
	(iMin,jMin)=np.nonzero(errorMatrix==np.min(errorMatrix))
	bestTheta = thetaVals[iMin]
	bestNoiseSigma = noiseSigmaVals[jMin]
	
	# Plot best-fitting curves:
	sliceDict['XVar'] = 'COn'
	sliceDict['theta'] = bestTheta
	sliceDict['noiseSigma'] = bestNoiseSigma
	plot1D(copy.copy(sliceDict), 'FC',saveResultDir =saveResultDir, whichRun =whichRun, tDel = 2000, tPen = 0, tND = tND, quickName = quickName,newFigure = 1)
	plot1D(copy.copy(sliceDict), 'RT',saveResultDir =saveResultDir, whichRun =whichRun, tDel = 2000, tPen = 0, tND = tND, quickName = quickName,newFigure = 1)
	
	if saveFig != 0:
		pl.figure(1)
		pl.savefig('/Users/Nick/Desktop/' + uniqueString + '_contour.eps')

	return bestTheta, bestNoiseSigma, errorMatrix



################################################################################
# Define the functions for the gaussian integral:
def intErfAB(a, b, mu=0, sigma=1):

	# Import necessary functions:
	import math
	
	# Ill need erf(x)
	def erf(x):
		# save the sign of x
		sign = 1
		if x < 0: 
			sign = -1
		x = abs(x)

		# constants
		a1 =  0.254829592
		a2 = -0.284496736
		a3 =  1.421413741
		a4 = -1.453152027
		a5 =  1.061405429
		p  =  0.3275911

		# A&S formula 7.1.26
		t = 1.0/(1.0 + p*x)
		y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
		return sign*y
	
	a=a*1.
	b=b*1.
	mu=mu*1.
	sigma=sigma*1.
	C1 = (b-mu)/(math.sqrt(2)*sigma)
	C2 = (a-mu)/(math.sqrt(2)*sigma)
	P1 = erf(C1)
	P2 = erf(C2)
	
	return .5*(P1-P2)

################################################################################
# Report the threshold ratio of a given slice:
def threshRatio(sliceDict,saveResultDir = './savedResults', whichRun = 0, quickName = -1,tND=350):

	# Import necessary stuff:
	from copy import copy
	import numpy as np



	# Gather relevent data:
	xVals, FCData = export1D(copy(sliceDict),'FC',saveResultDir=saveResultDir,whichRun=whichRun,tND=tND,quickName=quickName)
	xVals, RTData = export1D(copy(sliceDict),'RT',saveResultDir=saveResultDir,whichRun=whichRun,tND=tND,quickName=quickName)
	
	# Target half-ways in the range of psych/chrono functions:
	psyTar = .75
	chrTar = (max(RTData) + min(RTData))/2.0
	
	# Set up interpolating functions:
	def FCInt(x): return np.interp(x, xVals, FCData)
	def RTInt(x): return np.interp(x, xVals, RTData)
	
	# Set up functions for root finding:
	def FCIntRoot(x): return FCInt(x) - psyTar
	def RTIntRoot(x): return RTInt(x) - chrTar
	
	# Quick and dirty bisection method:
	def bisection(f,lb,ub):

		# Settings:
		eps = .001
		lb = float(lb)
		ub = float(ub)

		# Bisection method:
		while (abs(ub - lb) > 2*eps):
			mid = (ub + lb)/2
			if ((f(lb)*f(mid)) < 0):
				ub = mid
			elif ((f(ub)*f(mid)) < 0):
				lb = mid
			else:
				break
		return mid
	
	# Compute halfway-time-threshold (HTT) and halfway-accuracy-threshold (HAT)
	lb = min(xVals)
	ub = max(xVals)
	HAT = bisection(FCIntRoot,lb,ub)
	HTT = bisection(RTIntRoot,lb,ub)

	return HTT, HAT, HTT/HAT

################################################################################
# Plot 1-d Log plot:
def plot1DLog(sliceDict, whatToPlot, saveResultDir = './savedResults', whichRun = 0, quickName = -1, tDel = 2000, tPen = 0, tND=350, newFigure = 1, showFig = 1):

	# Import necessary stuff:
	import pylab as pl
	from copy import copy

	if whatToPlot == 'Both':
		plot1DLog(copy(sliceDict),'FC',saveResultDir=saveResultDir,whichRun=whichRun,quickName=quickName,tDel=tDel,tPen=tPen,tND=tND,newFigure = 1)
		plot1DLog(copy(sliceDict),'RT',saveResultDir=saveResultDir,whichRun=whichRun,quickName=quickName,tDel=tDel,tPen=tPen,tND=tND,newFigure = 1)
		
		return
	else:

		# Get data:
		X, Y = export1D(sliceDict, whatToPlot, tND=tND, whichRun=whichRun, quickName=quickName,saveResultDir=saveResultDir)
		
		# Determine if FD or not:
		inputSet = getSettings(saveResultDir=saveResultDir,whichRun=whichRun,quickName=quickName)[1]
		if inputSet == 0:
			inputSet = 'FR'
		else:
			inputSet = 'FD'
		
		# Get Psychophysics Data:
		roitmanFC, roitmanRT = getRoitmanPsyChr(inputSet)
		roitmanXPrime = pl.log(pl.array([3.2,6.4,12.8,25.6,51.2]))
		roitmanX = pl.concatenate(([0],roitmanXPrime - (roitmanXPrime[0]-(roitmanXPrime[-1]-roitmanXPrime[-2]))))

		# Make plotting data:
		logXPrime = pl.log(X[1:])
		logX = pl.concatenate(([0], logXPrime - (roitmanXPrime[0]-(roitmanXPrime[-1]-roitmanXPrime[-2]))))
		if whatToPlot == 'RT':
			roitmanY = roitmanRT
		elif whatToPlot == 'FC':
			roitmanY = roitmanFC

		# Plot data, with pooling noise and roitman dots:
		if showFig == 1:
			if newFigure: pl.figure()
			pl.plot(logX,Y)

			pl.plot(roitmanX, roitmanY,'o')

			# Set Axes:
			if whatToPlot == 'FC':
				XTick = pl.concatenate((pl.array([0]),pl.array([3.2,6.4,12.8,25.6,51.2])))
				L = min(logX)
				R = max(logX)
				M = 1 + 0.05*(1-.5)
				m = .5 - 0.05*(1-.5)
				pl.xlim(L,R)
				pl.ylim(m,M)
				pl.xticks(roitmanX, XTick)
				pl.yticks(pl.linspace(.5,1,6))
				pl.xlabel('Dot Coherence (C)')
				pl.ylabel('FC')
			elif whatToPlot == 'RT':
				XTick = pl.concatenate((pl.array([0]),pl.array([3.2,6.4,12.8,25.6,51.2])))
				YTicks=[400,600,800,1000]
				L = min(logX)
				R = max(logX)
				m=min(YTicks)
				M=max(YTicks)
				pl.xlim(L,R)
				pl.ylim(m,M)
				pl.xticks(roitmanX, XTick)
				pl.yticks(pl.floor(pl.linspace(min(YTicks),max(YTicks),len(YTicks))))
				pl.xlabel('Dot Coherence (C)')
				pl.ylabel('RT (ms.)')


	return ((roitmanX, roitmanY),(logX, Y))





