from __future__ import division		# allows '//' to mean integer division
import numpy
import Watershed as ws
#import scipy.special		# for the .kn() "modified Bessel function"
import scipy.signal		# for sepfir2d() and gaussian()


class ScaleSpace_Blob :
	def __init__(self) :
		self.isBright = None
		self.signifVal = None
		self.appearance = None
		self.disappearance = None

		# This will list grey blobs for scale each level.
		self.grey_blobs = []
		# This will list the support region for each scale level
		self.support_regions = []
		# This will list the scale value for each scale level
		self.scale_values = []
		
		self.approp_scalelvl = None
		self.approp_greyblob = None
		self.atBoundary = None

	def Add_Greylevel_Blobs(self, blobs, scaleVal) :
		new_support = ws.Support_Region()

		for aBlob in blobs :
			new_support.AddSupport(aBlob.support_region.pixels)

		self.support_regions.append(new_support)
		self.grey_blobs.append(blobs)
		self.scale_values.append(scaleVal)

	def Add_Greylevel_Blob(self, aBlob) :
		# Assumes that we are adding one blob to the current level.
		# Therefore, a level must exist first...
		self.support_regions[-1].AddSupport(aBlob.support_region.pixels)
		self.grey_blobs[-1].append(aBlob)


class Bifurcation_Event :
	def __init__(self, event_type = None, scaleblobs_above = [], scaleblobs_below = [], position = None, scaleVal = None) :
		self.event_type = event_type
		self.scaleblobs_above = scaleblobs_above
		self.scaleblobs_below = scaleblobs_below
		self.position = position
		self.scaleVal = scaleVal


class Scale_Level :
	def __init__(self, greyBlobs, image, scaleVal, coarserLevel = None, finerLevel = None) :
		self.scaleVal = scaleVal
		self.image = image
		self.greyBlobs = greyBlobs
#		self.greyBlobs_dark = []
#		self.coarserLevel = coarserLevel
#		self.finerLevel = finerLevel



class Primal_Sketch :
	def __init__(self) :
		self.scale_levels = []

		
		self.scaleBlobs_bright = []
		self.bifurcation_bright = []

		self.scaleBlobs_dark = []
		self.bifurcation_dark = []

	def CreateSketch(self, image, scale_levels) :
		if (len(self.scale_levels) == 0) :
			# No sketch exists, so start with the base 
			print "Base scale"
			greyblobs, basinMarks = ws.Watershed_Transform(image)
			self.AddBlobs(greyblobs, basinMarks, 0.0)
		

		for aScale in scale_levels :
			newImage = self.DoConvolve(image, aScale, (4 * (aScale // 2)) + 3).astype(int)
			print "At level: ", aScale, "  Image max:", newImage.max(), "   Image min:", newImage.min()
			greyblobs, basinMarks = ws.Watershed_Transform(newImage)
			self.AddBlobs(greyblobs, basinMarks, aScale)




	def DoConvolve(self, image, scale_level, winSize) :
		# NOTE: I know this isn't technically the best approach.
		# There is supposedly a better way using Bessel functions,
		#  but until I get a better idea how this is supposed to be
		#  implemented, I will do it this way.
		kernel = scipy.signal.gaussian(winSize, scale_level)
		return scipy.signal.sepfir2d(image, kernel / kernel.sum(),
						    kernel / kernel.sum())

	def AddBlobs(self, greyblobs, image, scaleVal) :
		# Right now, assume that the scale size is increasing monotonically.
		self.scale_levels.append(Scale_Level(greyblobs, image, scaleVal))

	def Create_ScaleSpace_Blobs(self) :
		# Skip out early if there is nothing to process
		if len(self.scale_levels) == 0 : return

		# Contains a list of the active scale space blobs
		# This gets modified at each scale level
		active_scale_blobs = []

		# We are going to loop through the levels (starting at the level
		#    that will likely have the fewest greylevel blobs).
		# We will loop through these levels with "atLevel".
		for atLevel in reversed(self.scale_levels) :
			
			# Map of grey blobs to scale blobs
			greyMap = [None] * len(atLevel.greyBlobs)
			# Map of scale blobs to grey blobs
			scaleMap = [None] * len(active_scale_blobs)

			for greyIndex, aGreyBlob in enumerate(atLevel.greyBlobs) :
				for blobIndex, aScaleBlob in enumerate(active_scale_blobs) :
					# The scalespace blob should be created with at least one level of data.
					# so accessing the last element of the support_regions list should be ok.
					sharedPixs = aGreyBlob.support_region.pixels.intersection(aScaleBlob.support_regions[-1].pixels)

					if len(sharedPixs) > 0 : 
						greyMap[greyIndex] = blobIndex
						scaleMap[blobIndex] = greyIndex


			# So, whichever entries in greyMap has "None", then we have a creation event
			# So, whichever entries in scaleMap has "None", then we have an extinction event
		
			# For each unique entry in greyMap we have a continuation event
			#     This should correspond to some unique entries in scaleMap.

			# For non-unique entries in greyMap, we have a scale blob associated with multiple grey blobs
			# For non-unique entries in scaleMap, we have a grey blob associated with multiple scale blobs

			# TODO: Much more work is needed, including the addition and utilization of bifurcation events
			#       Also need to include the merging of scale space blobs

			newLevelAdded = [False] * len(active_scale_blobs)
			for greyIndex, scaleBlobIndex in enumerate(greyMap) :
				if scaleBlobIndex is not None :
					if newLevelAdded[scaleBlobIndex] :
						active_scale_blobs[scaleBlobIndex].Add_Greylevel_Blob(atLevel.greyBlobs[greyIndex])
					else :
						active_scale_blobs[scaleBlobIndex].Add_Greylevel_Blobs([atLevel.greyBlobs[greyIndex]], atLevel.scaleVal)
						newLevelAdded[scaleBlobIndex] = True
					

			# NOTE: These two actions are done last because they modify the active_scale_blobs list
			for scaleBlobIndex in range(len(scaleMap) - 1, -1, -1) :
				if scaleMap[scaleBlobIndex] is None :
					# We have an extinction event!  Remove from the active list
					# TODO: Probably should modify the scaleblob a bit before eliminating it...
					del active_scale_blobs[scaleBlobIndex]


			for greyIndex, scaleBlobIndex in enumerate(greyMap) :
				if scaleBlobIndex is None :
					# We have a creation of a new scale space blob!
					new_blob = ScaleSpace_Blob()
					new_blob.Add_Greylevel_Blobs([atLevel.greyBlobs[greyIndex]], atLevel.scaleVal)
					active_scale_blobs.append(new_blob)
					self.scaleBlobs_bright.append(new_blob)
			

			



