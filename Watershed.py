
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
		for p in pixels[q] :
			if marked[p] == UNMARKED :
				isCenter = True
				#for aCenter in centers[q] :
				#	if aCenter.isNeighbor(p) :
				#		aCenter.addPoint(p)
				#		marked[p] = 1
	#			#		print q, "New Neighbor!"
				#		isCenter = False
				#		break
				
				markedSoFar = []
				for ii in range(max(p[0] - 1, 0), min(p[0] + 2, image.shape[0])) :
					for jj in range(max(p[1] - 1, 0), min(p[1] + 2, image.shape[1])) :
						if marked[ii, jj] == UNMARKED :
							# Only consider points that are lower than you right now.
							# Other points that are at the same level will be dealt with later.
							#if image[ii, jj] < q :
								marked[ii, jj] = 1
								markedSoFar.append((ii, jj))
						else :
							isCenter = False
							break

					if ~isCenter : break

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

	# for now... will later make this a parameter.
	#maxDepth = image.max()
	for delta in range(maxDepth) :
		deferredToNext = []
		print len(globs)
		for level in range(image.max(), -1, -1) :
			# Hysteresis level
			hlevel = level - delta
			print level, len(centers[level]), len(deferredToNext)

			centers[level] += deferredToNext
			deferredToNext = []

			foothills = []

			for centIndex, aCenter in enumerate(centers[level]) :
				if basins[aCenter.positions[0]] == UNMARKED :
					(basin, foothill) = Capture(basins, aCenter, basinNumber, image, hlevel, saliency)
					if basin is not None :
						foothills.append((foothill, centers[level][centIndex]))
						globs.append(basin)
						basinNumber += 1
					else :
						# Defer to next iteration to see if it will get big enough
						centers[level][centIndex].grey_val -= 1
						#deferredToNext.append(centers[level][centIndex])

			RemoveFoothills(basins, image, centers, foothills, hlevel)

	return globs, basins


def RemoveFoothills(basins, image, centers, foothills, hlevel) :
	for (foothill, center) in foothills :
		while len(foothill) > 0 :
			pixel = foothill.pop()
			basins[pixel] = GLOBBED

			# Checking the neighbors
			for ii in range(max(pixel[0] - 1, 0), min(pixel[0] + 2, image.shape[0])) :
				for jj in range(max(pixel[1] - 1, 0), min(pixel[0] + 2, image.shape[1])) :
					if basins[ii, jj] == UNMARKED :
						if ( image[ii, jj] >= 0 and image[ii, jj] < hlevel
						     and (image[ii, jj] <= image[pixel] or IsClosest(pixel, center, centers, image))) :
							foothill.append((ii, jj))



def IsClosest(pixel, center, centers, image) :
	mydist = min([(pixel[0] - centPos[0])**2 + (pixel[1] - centPos[1])**2
					for centPos in center.positions])

	binthresh = center.grey_val / 2

	for otherBin in range(binthresh, len(centers)) :
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
		for ii in range(max(pixel[0] - 1, 0), min(pixel[0] + 2, image.shape[0])) :
			for jj in range(max(pixel[1] - 1, 0), min(pixel[1] + 2, image.shape[1])) :
				if (basins[ii, jj] == UNMARKED) :
					if (~willBeConsideredAgain 
					    and image[ii, jj] >= 0 and image[ii, jj] < center.grey_val) :
						willBeConsideredAgain = True

					if image[ii, jj] >= hlevel :
						neighbors.append((ii, jj))
					elif image[ii, jj] >= 0 :
						# Not quite large enough to be a possible neighbor
						foothill.append((ii, jj))
	

	if len(basin.support_region.pixels) < saliency and willBeConsideredAgain :
		# Basin has not been captured
		# Now I need to undo what I have done...
		for p in basin.support_region.pixels :
			basins[p] = UNMARKED
		basin = None
		foothill = None

	return (basin, foothill)



