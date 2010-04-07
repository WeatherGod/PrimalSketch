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
		
		self.approp_scalelvl = None
		self.approp_greyblob = None
		self.atBoundary = None

	def Add_Greylevel_Blobs(self, blobs) :
		new_support = ws.Support_Region()

		for aBlob in blobs :
			new_support.AddSupport(aBlob.support_region.pixels)

		self.support_regions.append(new_support)
		self.grey_blobs.append(blobs)


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

	def Add_GreyLevel_Blobs(self, greyblobs, image, scaleVal) :
		# Right now, assume that the scale size is increasing monotonically.
		self.scale_levels.append(Scale_Level(greyblobs, image, scaleVal))

		# Map of greys to scales
		greyMap = [None] * len(greyblobs)
		# Map of scales to greys
		scaleMap = [None] * len(self.scaleBlobs_bright)

		# Oops, I need to make sure I only grab scale blobs that are in the previous level,
		#	I should skip all of the extinct scale blobs.
		for greyIndex, aGreyBlob in enumerate(greyblobs) :
			for blobIndex, aScaleBlob in enumerate(self.scaleBlobs_bright) :
				# The scalespace blob should be created with at least one level of data.
				# so accessing the last element of the support_regions list should be ok.
				sharedPixs = aGreyBlob.support_region.pixels.intersection(aScaleBlob.support_regions[-1].pixels)

				if len(sharedPixs) > 0 : 
					greyMap[greyIndex] = blobIndex
					scaleMap[blobIndex] = greyIndex

		# So, whichever entries in greyMap has "None", then we have a creation event
		#     whichever entries in scaleMap has "None", then we have an extinction event

		# For each unique entry in greyMap we have a continuation event
		#     This should correspond to some unique entries in scaleMap.

		# For non-unique entries in greyMap, we have a scale blob associated with multiple grey blobs
		# For non-unique entries in scaleMap, we have a grey blob associated with multiple scale blobs

			
			



