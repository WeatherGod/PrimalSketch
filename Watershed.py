
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
		self.positions = pos
		self.grey_val = val
		self.grey_blob = greyblob

	def isNeighbor(self, pos) :
		return(min([(pos[0] - aPoint[0])**2 + (pos[1] - aPoint[1])**2
				for aPoint in self.positions]) <= 1.0)

	def addPoint(self, pos) :
		self.positions.append(pos)

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


def Neighbors(pos, rangeVal, shape) :
	neighbors = []
	for x in range(max(pos[0] - rangeVal, 0), min(pos[0] + rangeVal + 1, shape[0])) :
		for y in range(max(pos[1] - rangeVal, 0), min(pos[1] + rangeVal + 1, shape[1])) :
			if (x != pos[0] or y != pos[1]) :
				neighbors.append((x, y))

	return neighbors


def Watershed_Transform(image, maxDepth, saliency) :
	# Assume that image is a normalized, quantized image...
	pixels = [ [] for i in range(image.max() + 1)]

	# loading the 'pixels' array with coordinates of pixels with associated values
	for x in range(image.shape[0]) :
		for y in range(image.shape[1]) :
			pixels[image[x, y]].append((x, y))

	centers = [ [] for i in range(image.max() + 1)]
	marked = UNMARKED * numpy.ones_like(image)
	for q in range(image.max(), -1, -1) :
		"""
		for p in pixels[q] :
			if marked[p] == UNMARKED :
				isCenter = True
				for aCenter in centers[q] :
					if aCenter.isNeighbor(p) :
						aCenter.addPoint(p)
	#			#		print q, "New Neighbor!"
						isCenter = False
						break

				if isCenter :
					centers[q].append(Extremum_Region([p], q))

				marked[p] = 1
		"""
		for p in pixels[q] :
			if marked[p] == UNMARKED :
				isCenter = True
				markedSoFar = []
				for point in Neighbors(p, 1, image.shape) :
					if marked[point] == UNMARKED :
						marked[point] = 1
						markedSoFar.append(point)
					else :
						# p touches an already marked point,
						# so it can't be a center, but also don't
						# want to bother with it again.
						marked[p] = 1
						isCenter = False
						break


				if isCenter : 
					centers[q].append(Extremum_Region([p], q))
					marked[p] = 1
#					print q, "New Center!"
				else :
					# time to undo the markings
					for aPos in markedSoFar :
						marked[aPos] = UNMARKED

#		print q, len(centers[q])

	
	# Zero for background, -1 for unchecked, positive values for blob number
	basins = UNMARKED * numpy.ones_like(image)
	# Initializing the basin number
	basinNumber = 1

	
	globs = []

	for delta in range(maxDepth) :
		deferredToNext = []
		print len(globs)
		for level in range(image.max(), -1, -1) :
			# Hysteresis level
			hlevel = level - delta
			#print level, len(centers[level]), len(deferredToNext)

			centersTmp = centers[level] + deferredToNext
			deferredToNext = []

			foothills = []

			for centIndex, aCenter in enumerate(centersTmp) :
				if basins[aCenter.positions[0]] == UNMARKED :
					(basin, foothill) = Capture(basins, aCenter, basinNumber, image, hlevel, saliency)
					if basin is not None :
						foothills.append((foothill, centersTmp[centIndex]))
						globs.append(basin)
						basinNumber += 1
					else :
						# Defer to next iteration to see if it will get big enough
						centersTmp[centIndex].grey_val -= 1
						deferredToNext.append(centersTmp[centIndex])

			print "%3d Centers: %4d Deferred: %3d  Globs: %4d  Foothills: %4d" %  (level, len(centers[level]), len(deferredToNext), len(globs), len(foothills))
			RemoveFoothills(basins, image, centers, foothills, hlevel)


	return globs, basins


def RemoveFoothills(basins, image, centers, foothills, hlevel) :
	for (foothill, center) in foothills :
		while len(foothill) > 0 :
			pixel = foothill.pop()
			basins[pixel] = GLOBBED

			# Checking the neighbors
			for point in Neighbors(pixel, 1, image.shape) :
				if basins[point] == UNMARKED :
					if ( image[point] >= 0 and image[point] < hlevel
					     and (image[point] <= image[pixel] or IsClosest(pixel, center, centers, image))) :
						foothill.append(point)



def IsClosest(pixel, center, centers, image) :
	mydist = min([(pixel[0] - centPos[0])**2 + (pixel[1] - centPos[1])**2
					for centPos in center.positions])

	binthresh = center.grey_val / 2

	for otherBin in range(binthresh, len(centers)) :
		#print len(centers[otherBin])
		for otherCenter in centers[otherBin] :
			if mydist > min([(pixel[0] - centPos[0])**2 + (pixel[1] - centPos[1])**2 for centPos in otherCenter.positions]) :
				return False

	return True


def Capture(basins, center, basinNumber, image, hlevel, saliency) :
	neighbors = []

	foothill = []
	basin = GreyLevel_Blob()
	basin.extremum = center
	basin.support_region = Support_Region()
	basin.support_region.pixels = []

	willBeConsideredAgain = False
	
	neighbors += center.positions
	while (len(neighbors) > 0) :
		pixel = neighbors.pop()
		if basins[pixel] != UNMARKED : continue	# already processed

		basins[pixel] = basinNumber
		basin.support_region.pixels.append(pixel)

		# Checking the neighbors
		for point in Neighbors(pixel, 1, image.shape) :
			if (basins[point] == UNMARKED) :
				if (~willBeConsideredAgain 
				    and image[point] >= 0 and image[point] < center.grey_val) :
					willBeConsideredAgain = True

				if image[point] >= hlevel :
					neighbors.append(point)
				elif image[point] >= 0 :
					# Not quite large enough to be a possible neighbor
					foothill.append(point)
	

	if len(basin.support_region.pixels) < saliency and willBeConsideredAgain :
		# Basin has not been captured
		# Now I need to undo what I have done...
		for p in basin.support_region.pixels :
			basins[p] = UNMARKED
		basin = None
		foothill = None

	return (basin, foothill)



