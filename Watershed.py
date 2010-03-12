
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

GLOBBED = -3
UNMARKED = -1


def Neighbors(pos, rangeVal, shape, inclSelf = False) :
	neighbors = []
	for x in range(max(pos[0] - rangeVal, 0), 
		       min(pos[0] + rangeVal + 1, shape[0])) :
		for y in range(max(pos[1] - rangeVal, 0), 
			       min(pos[1] + rangeVal + 1, shape[1])) :
			if (inclSelf or x != pos[0] or y != pos[1]) :
				neighbors.append((x, y))

	return neighbors


def FindCenters(image) :
	# Assume that image is a normalized, quantized image...
	pixels = [ [] for i in range(image.max() + 1)]

	# loading the 'pixels' array with coordinates of pixels with associated values
	for x in range(image.shape[0]) :
		for y in range(image.shape[1]) :
			pixels[image[x, y]].append((x, y))

	return pixels
	


def Watershed_Transform(image, maxDepth) :


	centers = FindCenters(image)

	
	# Zero for background, -1 for unchecked, positive values for blob number
	basins = UNMARKED * numpy.ones_like(image)
	# Initializing the basin number
	basinNumber = 0

	
	globs = []

	deferredToNext = []
	#print delta, len(globs)
	for level in range(image.max(), -1, -1) :
		# Hysterisis level.  Don't let it get below 0.
		#hlevel = max(level - maxDepth, 0)
		#print level, len(centers[level]), len(deferredToNext)

		#centersTmp = centers[level] + deferredToNext
		#deferredToNext = []

		markedSoFar = {}
		clearOut = []

		for aPix in centers[level] :
			connectedTo = set([basins[neighPix] for neighPix in Neighbors(aPix, 1, image.shape)])
			connectedTo -= set([UNMARKED, GLOBBED])
				
			if len(connectedTo) == 0 :
				# Doesn't touch any existing basins, must be a new basin!
				basins[aPix] = basinNumber
				markedSoFar[basinNumber] = [aPix]
				basinNumber += 1
			elif len(connectedTo) >= 1 :
				# only touches one known basin, so let that basin grow!
				thisBasin = connectedTo.pop()
				basins[aPix] = thisBasin
				markedSoFar[thisBasin].append(aPix)
				
				if len(connectedTo) > 1 :
					# Ah, this isopleth is the limit for the connectedTo basins!
					clearOut += connectedTo

		for basinToClean in clearOut :
			for aPix in markedSoFar[basinToClean] :
				basins[aPix] = GLOBBED
				

#				(basin, captured) = Capture(image, basins, aCenter, basinNumber, hlevel)
#				if not captured :
#					# Defer to next iteration to see if it will get big enough
#					centersTmp[centIndex].grey_val -= 1
#					deferredToNext.append(centersTmp[centIndex])
#				elif basin is not None :
#					globs.append(basin)
#					basinNumber += 1
#
#			
#		#print "%3d  Centers: %4d  Deferred: %3d  Globs: %4d  Foothills: %4d  " %  (level, len(centers[level]), len(deferredToNext), len(globs), len(foothills))
#		RemoveFoothills(image, basins, hlevel, centers, foothills)


	return globs, basins


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



