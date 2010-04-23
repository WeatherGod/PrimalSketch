from __future__ import division		# allows '//' to mean integer division
import numpy
import Watershed as ws
#import scipy.special		# for the .kn() "modified Bessel function"
import scipy.signal		# for sepfir2d() and gaussian()


class ScaleSpace_Blob :
	def __init__(self, scale_level, idNum) :
		self.isBright = None
		self.signifVal = None
		self.appearance = scale_level
		self.disappearance = None
		self.id_num = idNum

		# This will list grey blobs for scale each level.
		self.grey_blobs = [ [] ]
		# This will list the support region for each scale level
		self.support_regions = [ws.Support_Region()]
		# This will list the scale value for each scale level
		self.scale_levels = [scale_level]
		
		self.approp_scalelvl = None
		self.approp_greyblob = None
		self.atBoundary = None

	def Add_Greylevel_Blobs(self, blobs, scale_lev) :
		new_support = ws.Support_Region()

		for aBlob in blobs :
			aBlob.scalespace_blob = self
			new_support.AddSupport(aBlob.support_region.pixels)

		self.support_regions.append(new_support)
		self.grey_blobs.append(blobs)
		self.scale_levels.append(scale_lev)

	def Add_Greylevel_Blob(self, aBlob) :
		# Assumes that we are adding one blob to the current level.
		# Therefore, a level must exist first...
		aBlob.scalespace_blob = self
		self.support_regions[-1].AddSupport(aBlob.support_region.pixels)
		self.grey_blobs[-1].append(aBlob)




class Bifurcation_Event :
	def __init__(self, event_type = None, scaleblobs_above = [], scaleblobs_below = [], position = None, scale_lev = None) :
		self.event_type = event_type
		self.scaleblobs_above = scaleblobs_above
		self.scaleblobs_below = scaleblobs_below
		self.position = position
		self.scale_level = scale_lev


class Scale_Level :
	def __init__(self, greyBlobs, image, scaleVal, coarserLevel = None, finerLevel = None) :
		self.scaleVal = scaleVal
		self.image = image
		self.greyBlobs = greyBlobs
#		self.greyBlobs_dark = []
#		self.coarserLevel = coarserLevel
#		self.finerLevel = finerLevel
		#self.greyMap = [None] * len(greyBlobs)



class Primal_Sketch :
	def __init__(self) :
		self.scale_levels = []

		self.scaleBlobs_bright = []
		self.bifurcation_bright = []

		self.scaleBlobs_dark = []
		self.bifurcation_dark = []
		#self.grey2scale_maps = []

		self.scaleBlob_Marks = []
		self.UNMARKED = -1



	def CreateSketch(self, image, scale_values) :
		#if (len(self.scale_levels) == 0) :
		#	# No sketch exists, so start with the base 
		#	print "Base scale"
		#
		#	greyblobs, basinMarks = ws.Watershed_Transform(image)
		#	self.AddBlobs(greyblobs, basinMarks, 0.0)
		

		for aScale in scale_values :
			if aScale == 0 :
				newImage = image
			else :
				newImage = self.DoConvolve(image, aScale, (4 * (aScale // 2)) + 3).astype(int)

			print "At level: ", aScale #, "  Image max:", newImage.max(), "   Image min:", newImage.min()
			
			greyblobs, basinMarks = ws.Watershed_Transform(newImage)
			self.AddBlobs(greyblobs, basinMarks, aScale)
			self.Create_ScaleSpace_Blobs()




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
		self.scaleBlob_Marks.append( numpy.zeros(image.shape) + self.UNMARKED )
		

	def Create_ScaleSpace_Blobs(self) :
		# Skip out early if there is nothing to process
		if len(self.scale_levels) == 0 : return

		currIDNum = len(self.scaleBlobs_bright)

		atLevel = self.scale_levels[-1]
		
		if len(self.scale_levels) > 1 :
			previousBlobMarks = self.scaleBlob_Marks[-2]
		else :
			previousBlobMarks = numpy.zeros(atLevel.image.shape) + self.UNMARKED
		
		ignoreThese = frozenset([self.UNMARKED])
		scaleBlobsMatched = set([])

		for greyIndex, aGreyBlob in enumerate(atLevel.greyBlobs) :
			blobIndices = set([previousBlobMarks[anIndex] for anIndex in aGreyBlob.support_region.pixels]) - ignoreThese

			# So, if len(blobIndices) is zero, then we have a new scalespace blob!
			#     if len(blobIndices) is one, then we have a continuation,
			#     if len(blobIndices) is greater than one, we have a bifurcation event.

			if len(blobIndices) == 0 :

				for anIndex in aGreyBlob.support_region.pixels :
					self.scaleBlob_Marks[-1][anIndex] = currIDNum

				new_blob = ScaleSpace_Blob(atLevel, currIDNum)
				currIDNum += 1
				new_blob.Add_Greylevel_Blob(aGreyBlob)
				self.scaleBlobs_bright.append(new_blob)

			elif len(blobIndices) == 1 :
				blobIndex = list(blobIndices)[0]
				print blobIndex
				for anIndex in aGreyBlob.support_region.pixels :
					self.scaleBlob_Marks[-1][anIndex] = blobIndex

				if blobIndex in scaleBlobsMatched :
					self.scaleBlobs_bright[blobIndex].Add_GreyLevel_Blob(aGreyBlob)
				else :
					self.scaleBlobs_bright[blobIndex].Add_GreyLevel_Blobs([aGreyBlob], atLevel)
					scaleBlobsMatched.update([blobIndex])

			else :
				# TODO: Figure out what to do here...
				print "Grrr... a bad component labeling issue.  Just pick the first one and go with it."
				blobIndex = list(blobIndices)[0]
				for anIndex in aGreyBlob.support_region.pixels :
					self.scaleBlob_Marks[-1][anIndex] = blobIndex

                                if blobIndex in scaleBlobsMatched :
                                        self.scaleBlobs_bright[blobIndex].Add_GreyLevel_Blob(aGreyBlob)
                                else :
                                        self.scaleBlobs_bright[blobIndex].Add_GreyLevel_Blobs([aGreyBlob], atLevel)
                                        scaleBlobsMatched.update([blobIndex])


		# Any scale blobs from the previous run that are unmatched needs to be discontinued.
		unmatchedBlobs = set(previousBlobMarks.flat) - scaleBlobsMatched - ignoreThese
		for blobIndex in unmatchedBlobs :
			self.scaleBlobs_bright[blobIndex].disappearance = atLevel
			
