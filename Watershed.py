
import numpy

class GreyLevel_Blob :
	def __init__(self) :
		self.polarity = True
		self.scale_level = None
		self.extremum = None
		self.saddle = None
		self.support_region = None



class Extremum_Region :
	def __init__(self, pos = [], val = None, greyblob = None) :
		self.position = pos
		self.grey_val = val
		self.grey_blob = greyblob

#	def isNeighbor(self, pos) :
#		return(min([(pos[0] - aPoint[0])**2 + (pos[1] - aPoint[1])**2
#				for aPoint in self.positions]) <= 1.0)
#
#	def addPoint(self, pos) :
#		self.positions.append(pos)

class Saddle_Region :
	def __init__(self, pos = [], val = None, greyblobs = []) :
		self.position = pos
		self.grey_val = val
		self.grey_blobs = greyblobs


class Support_Region :
	def __init__(self) :
		self.blob_area = None
		self.first_moment = None
		self.second_moment = None
		self.pixels = None
		self.atBoundary = None


UNMARKED = -1
GLOBBED = -3

def Neighbors(pos, rangeVal, shape, inclSelf = False) :
	neighbors = []
	for x in range(max(pos[0] - rangeVal, 0), 
		       min(pos[0] + rangeVal + 1, shape[0])) :
		for y in range(max(pos[1] - rangeVal, 0), 
			       min(pos[1] + rangeVal + 1, shape[1])) :
			if (inclSelf or x != pos[0] or y != pos[1]) :
				neighbors.append((x, y))

	return neighbors


def SortImage(image) :
	# Assume that image is a normalized, quantized image...
	pixels = [ [] for i in range(image.max() + 1)]

	# loading the 'pixels' array with coordinates of pixels with associated values
	for x in range(image.shape[0]) :
		for y in range(image.shape[1]) :
			pixels[image[x, y]].append((x, y))

	pixels.reverse()
	return pixels

def ConnectedComponets(pixels, imageShape) :
	"""
	Process the pixels of the same gray level value.  Connect the neighboring pixels
	"""
	components = []
	componentID = 0

	componentMarks = UNMARKED * numpy.ones(imageShape)

	for aPix in pixels :
		theNeighbors = Neighbors(aPix, 1, imageShape)
		connectedIsos = set([componentMarks[neighPix] for neighPix in theNeighbors]) - set([UNMARKED])

		
		if len(connectedIsos) == 0 :
			# A lonely pixel!  Start a new isopleth!
			components.append([aPix])
			componentmarks[aPix] = componentID
			componentID += 1
		else :
			# Ah, this pixel neighbors at least one existing isopleth!
			thisIso = connectedIsos.pop()

			# If there are any other connected Isopleths, then we have
			# a degenerate signal issue!  Therefore, we are going to have
			# to merge these other isopleths into the first one.
			for otherIso in connectedIsos :
				for otherPix in components[otherIso] :
					componentMarks[otherPix] = thisIso
				components[thisIso] += components[otherIso]
				components[otherIso] = None

			componentMarks[aPix] = thisIso
			components[thisIso].append(aPix)

	# Squeeze out the empy isopleth entries
	components = [aPleth for aPleth in components if aPleth is not None]

	return components
	


def GetIsopleths(pixels, isoMarks, isoplethID) :
	print len(pixels), isoplethID
	curisoID = isoplethID
	isopleths = []

	componentMap = []

	
	# Processing all pixels for an isolevel
	for aPix in pixels :
		allComps = set([isoMarks[neighPix] for neighPix in Neighbors(aPix, 1, isoMarks.shape)]) - set([UNMARKED, GLOBBED])

		# Return only neighboring isopleth IDs of those isopleths in pixels
		connectedComponents = allComps - set(range(curisoID))

		# Return only neighboring isopleth IDs of higher isopleths
		connectedHighers = allComps - set(range(curisoID, isoplethID))

		if len(connectedComponents) == 0 :
			# A lonely pixel!  Start a new isopleth!
			isopleths.append({'higherIsos': set(connectedHighers), 
					  'pixels': [aPix],
					  'componentID': isoplethID,
					  'components': set([isoplethID])})
			isoMarks[aPix] = isoplethID
			componentMap.append(isoplethID)
			isoplethID += 1
		else :

			connectedIsos = set([componentMap[aComponent - curisoID] for aComponent in connectedComponents])

			# Ah, this pixel neighbors at least one existing isopleth!
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

	print "Merging..."
	# Merging the components, working backwards
	for index in range(len(isopleths) - 1, -1, -1) :
		if componentMap[index] != isopleths[index]['componentID'] :
			# This item has to get moved!
			isopleths[componentMap[index] - curisoID]['higherIsos'].update(isopleths[index]['higherIsos'])
			isopleths[componentMap[index] - curisoID]['pixels'] += isopleths[index]['pixels']
			isopleths[componentMap[index] - curisoID]['components'].update(isopleths[index]['components'])

			isopleths[index] = None

	print "Consolidating..."	
	finalIsopleths = [isopleths[anIsoID - curisoID] for anIsoID in set(componentMap)]
	for index, anIso in enumerate(finalIsopleths) :
		anIso['componentID'] = index + curisoID
		for aPix in anIso['pixels'] :
			isoMarks[aPix] = index + curisoID

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
			globs.append({'componentIDs': [componentNum], 'start_level': level, 'stop_level': None})

			
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
					if globs[basinToAssign]['stop_level'] is not None :
						basinToAssign = GLOBBED
					else :
						globs[basinToAssign]['componentIDs'].append(componentNum)
				else :
					# Oops, we got multiple basins!  Gotta set them all to stop growing
					basinToAssign = GLOBBED

					
			# This is a catch-all to make sure that the globs get stopped if they have to be
			if basinToAssign == GLOBBED :
				for aBasin in basinsTouch :
					if aBasin != GLOBBED and globs[aBasin]['stop_level'] is None :
						globs[aBasin]['stop_level'] = level

		componentMap[componentNum] = basinToAssign
		for aPix in pixels :
			basinMarks[aPix] = basinToAssign

	return basinNumber

def ColorComponents(components, imageShape) :
	isoImage = UNMARKED * numpy.ones(imageShape, dtype = int)

	print components[0]['componentID'], components[-1]['componentID']

	for anIso in components :
		for aPix in anIso['pixels'] :
			isoImage[aPix] = anIso['componentID']

	return isoImage


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

	basins = []

	# We now process the isopleths, starting with the highest valued isopleths
	for index, components in enumerate(isopleths) :
		print index
		basinNumber = MarkComponents(components, basinMarks, basinNumber, len(pixels) - index, globs, componentMap)
		basins.append(basinMarks.copy())


	return globs, basins


"""

def RemoveFoothills(image, basins, hlevel, centers, foothills) :
	for (foothill, center) in foothills :
		while len(foothill) > 0 :
			pixel = foothill.pop()
			basins[pixel] = GLOBBED

			# Checking the neighbors
			for point in Neighbors(pixel, 1, image.shape) :
				if basins[point] == UNMARKED :
					if ( image[point] >= 0 and image[point] < hlevel
					     and (image[point] <= image[pixel] or IsClosest(pixel, center, centers))) :
						foothill.append(point)



def IsClosest(pixel, center, centers) :
	mydist = (pixel[0] - center.position[0])**2 + (pixel[1] - center.position[1])**2

	binthresh = center.grey_val / 2

	for otherBin in range(binthresh, len(centers)) :
		#print len(centers[otherBin])
		for otherCenter in centers[otherBin] :
			if mydist > ((pixel[0] - otherCenter.position[0])**2 + (pixel[1] - otherCenter.position[1])**2) :
				return False

	return True


def Capture(image, basins, center, basinNumber, hlevel, saliency, foothills) :
	neighbors = [center.position]

	foothill = []
	markedsofar = []

	willBeConsideredAgain = False
	
	while (len(neighbors) > 0) :
		pixel = neighbors.pop()
		if basins[pixel] != UNMARKED : continue	# already processed

		basins[pixel] = basinNumber
		markedsofar.append(pixel)

		# Checking the neighbors
		for point in Neighbors(pixel, 1, image.shape) :
			if (basins[point] == UNMARKED) :
				if (not willBeConsideredAgain 
				    and image[point] >= 0 and image[point] < center.grey_val) :
					willBeConsideredAgain = True

				if image[point] >= hlevel :
					neighbors.append(point)
				elif image[point] >= 0 :
					# Not quite large enough to be a possible neighbor
					foothill.append(point)

	if center.grey_val == 0 :
		willBeConsideredAgain = False

	bigEnough = (len(markedsofar) >= saliency)

	basin = None
	
	if bigEnough :
		foothills.append((foothill, center))
		basin = GreyLevel_Blob()
		basin.extremum = center
		basin.support_region = Support_Region()
		basin.support_region.pixels = markedsofar
	elif willBeConsideredAgain :
		# Basin has not been captured
		# Now I need to undo what I have done...
		for p in markedsofar :
			basins[p] = UNMARKED
		markedsofar = None
		foothill = None
		neighbors = None
	else :
		#print "Not being deferred!"
		# So, it is not big enough, and it won't be considered again,
		# Then mark them as globbed, so they won't cause confusion.
		for p in markedsofar :
			basins[p] = GLOBBED
		markedsofar = None
		foothill = None
		neighbors = None

	return (basin, (bigEnough or not willBeConsideredAgain))

"""

