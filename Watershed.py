
import numpy

class GreyLevel_Blob :
	def __init__(self, pixels, value, polarity = True) :
		self.polarity = polarity
		#self.scale_level = scaleLevel
		#self.scalespace_blob = None

		self.extremum = Extremum_Region(pixels, value, self)
		self.saddle = None
		self.support_region = Support_Region(pixels)

	def AddSupport(self, pixels) :
		self.support_region.AddSupport(pixels)





class Support_Region(set) :
	def __init__(self, pixels = []) :
#		self.first_moment = None
#		self.second_moment = None
		self.atBoundary = None
		set.__init__(self, pixels)

	def AddSupport(self, pixels) :
		self.update(pixels)

	def blob_area(self) :
		return len(self)

	def first_moment(self) :
		return list(numpy.average(list(self), axis=0))

	def second_moment(self) :
		return list(numpy.std(list(self), axis=0))

class Extremum_Region(Support_Region) :
	def __init__(self, pos = [], val = None, greyblob = None) :
		self.grey_val = val
		self.grey_blob = greyblob
		set.__init__(self, pos)

class Saddle_Region(Support_Region) :
	def __init__(self, pos = [], val = None, greyblobs = []) :
		self.grey_val = val
		self.grey_blobs = greyblobs
		set.__init__(self, pos)

UNMARKED = -1
GLOBBED = -3

def Neighbors(pos, rangeVal, shape, inclSelf = False) :
#	neighbors = []
	startx = pos[1] - rangeVal if (pos[1] - rangeVal) > 0 else 0
	starty = pos[0] - rangeVal if (pos[0] - rangeVal) > 0 else 0
	endx = pos[1] + rangeVal + 1 if (pos[1] + rangeVal + 1) < shape[1] else shape[1]
	endy = pos[0] + rangeVal + 1 if (pos[0] + rangeVal + 1) < shape[0] else shape[0]
	if not inclSelf :
		return [(y, x) for x in range(startx, endx) for y in range(starty, endy) if (y != pos[0] or x != pos[1])]
	else :
		return [(y, x) for x in range(startx, endx) for y in range(starty, endy)]


#	for x in range(startx, endx) :
#		for y in range(starty, endy) :
#			if (inclSelf or x != pos[0] or y != pos[1]) :
#				neighbors.append((x, y))

#	return neighbors


def SortImage(image) :
	# Assume that image is a normalized, quantized image...
	pixels = [ [] for i in range(image.max() + 1)]

	# loading the 'pixels' array with coordinates of pixels with associated values
	#for index, imageVal in numpy.ndenumerate(image) :
	#	pixels[imageVal].append(index)
	xIndices = range(image.shape[1])
	yIndices = range(image.shape[0])

	# Oddly enough, going through the indices is faster than using ndenumerate()
	for y in yIndices :
		for x in xIndices :
			pixels[image[y, x]].append((y, x))

	# Order from greatest pixel values to least.
	pixels.reverse()
	return pixels



def GetIsopleths(pixels, isoMarks, isoplethID) :
	#print len(pixels), isoplethID
	curisoID = isoplethID
	isopleths = []

	componentMap = []
	establishedComponents = frozenset(range(curisoID))
	nonComponents = frozenset([UNMARKED, GLOBBED])
	newComponents = set([])

	
	# Processing all pixels for an isolevel
	for aPix in pixels :
		allComps = set([isoMarks[neighPix] for neighPix in Neighbors(aPix, 1, isoMarks.shape)]) - nonComponents

		# Return only neighboring isopleth IDs of those isopleths in pixels
		connectedComponents = allComps - establishedComponents

		# Return only neighboring isopleth IDs of higher isopleths
		connectedHighers = allComps - newComponents

		if len(connectedComponents) == 0 :
			# A lonely pixel!  Start a new isopleth!
			isopleths.append({'higherIsos': connectedHighers.copy(), 
					  'pixels': [aPix],
					  'componentID': isoplethID,
					  'components': set([isoplethID])})
			isoMarks[aPix] = isoplethID
			componentMap.append(isoplethID)
			newComponents.update([isoplethID])
			isoplethID += 1
		else :
			# Use the component map to reduce this list.
			connectedIsos = set([componentMap[aComponent - curisoID] for aComponent in connectedComponents])

			# This pixel neighbors at least one existing isopleth.
			thisIso = min(connectedIsos)
			connectedIsos.remove(thisIso)

			# If there are any other connected Isopleths, then we have
			# a degenerate signal issue!  Therefore, we are going to have
			# to merge these other isopleths into the first one.
			if len(connectedIsos) > 0 :
				UpdateComponentMap(componentMap, connectedComponents, isopleths, curisoID, thisIso)
			
			isoMarks[aPix] = thisIso
			isopleths[thisIso - curisoID]['higherIsos'].update(connectedHighers)
			isopleths[thisIso - curisoID]['pixels'].append(aPix)
			isopleths[thisIso - curisoID]['components'].update(connectedComponents)

	print "Merging...", curisoID
	# Merging the components, working backwards
	for index in range(len(isopleths) - 1, -1, -1) :
		if componentMap[index] != isopleths[index]['componentID'] :
			# This item has to get moved!
			isopleths[componentMap[index] - curisoID]['higherIsos'].update(isopleths[index]['higherIsos'])
			isopleths[componentMap[index] - curisoID]['pixels'] += isopleths[index]['pixels']
			isopleths[componentMap[index] - curisoID]['components'].update(isopleths[index]['components'])

			isopleths[index] = None

	#print "Consolidating..."	
	finalIsopleths = [isopleths[anIsoID - curisoID] for anIsoID in set(componentMap)]
	for index, anIso in enumerate(finalIsopleths) :
		anIso['componentID'] = index + curisoID
		for aPix in anIso['pixels'] :
			isoMarks[aPix] = anIso['componentID']

	return finalIsopleths, curisoID + len(finalIsopleths)

def UpdateComponentMap(componentMap, componentList, isopleths, curisoID, thisIso) :
	for aComponent in componentList :
		if componentMap[aComponent - curisoID] != thisIso :
			componentMap[aComponent - curisoID] = thisIso
			UpdateComponentMap(componentMap, isopleths[aComponent - curisoID]['components'], 
					   isopleths, curisoID, thisIso)


def MarkComponents(components, basinMarks, basinNumber, level, globs, componentMap) :

	# Now, process each isopleth and figure out what it should be assigned as
	for anIso in components :
		# Ignore Globbed for the moment
		touchHigherIsos = anIso['higherIsos']
		pixels = anIso['pixels']
		componentNum = anIso['componentID']

		#print len(touchBasins), len(pixels)

		# By default, an isopleth will be called "Globbed"
		# Only if the isopleth is isolated or only touching
		# one basin will it get a basin number
		basinToAssign = GLOBBED

		if len(touchHigherIsos) == 0 :
			# A lonely isopleth!
			# Give it a new basin
			basinToAssign = basinNumber
			basinNumber += 1
			globs.append(GreyLevel_Blob(pixels, level))
			#print "Iso pixel len:",  len(pixels)
		else :
			# There is at least one higher isopleth touching this isopleth.
			# We need to know how many basins that represents...
			basinsTouch = set([componentMap[aComponent] for aComponent in touchHigherIsos])

			if GLOBBED in basinsTouch :
				# If even one of the touching basins is globbed,
				# then this one gets globbed.
				basinToAssign = GLOBBED
				
			else :
				if len(basinsTouch) == 1 :
					# Looking good, just need to check to see if the basin
					# is still allowed to grow.
					basinToAssign = basinsTouch.pop()
					if globs[basinToAssign].saddle is not None :
						# This blob has stopped growing
						basinToAssign = GLOBBED
					else :
						# This blob is still growing.
						globs[basinToAssign].AddSupport(pixels)
						
				else :
					# Oops, we got multiple basins!  Gotta set them all to stop growing
					basinToAssign = GLOBBED

					
			# This is a catch-all to make sure that the globs get stopped if they have to be
			if basinToAssign == GLOBBED :
				theGlobs = [globs[aBasin] for aBasin in basinsTouch if aBasin != GLOBBED 
											and globs[aBasin].saddle is None]
				newSaddle = Saddle_Region(pixels, level, theGlobs)
				for aBlob in theGlobs :
					aBlob.saddle = newSaddle

		componentMap[componentNum] = basinToAssign
		for aPix in pixels :
			basinMarks[aPix] = basinToAssign

	return basinNumber


def Watershed_Transform(image) :

	
	pixels = SortImage(image)

	isopleths = []
	isoMarks = UNMARKED * numpy.ones(image.shape, dtype=int)
	componentID = 0

	# Do all but the last level, since it will be background,
	# and typically has the most pixels of any isolevel.
	for level in range(len(pixels) - 1) :
		components, componentID = GetIsopleths(pixels[level], isoMarks, componentID)
		isopleths.append(components)
		
	globs = []
	componentMap = GLOBBED * numpy.ones(componentID + 1, dtype=int)

	# GLOBBED for background, UNMARKED for unchecked, positive values for blob/basin number
	basinMarks = UNMARKED * numpy.ones_like(image)
	# Initializing the basin number
	basinNumber = 0

	#basins = []

	# We now process the isopleths, starting with the highest valued isopleths
	for index, components in enumerate(isopleths) :
		#print index
		basinNumber = MarkComponents(components, basinMarks, basinNumber, len(pixels) - index, globs, componentMap)
		#basins.append(basinMarks.copy())

	return globs, basinMarks


