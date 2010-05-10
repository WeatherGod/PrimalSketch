from __future__ import division		# allows '//' to mean integer division
import numpy
import Watershed as ws
#import scipy.special		# for the .kn() "modified Bessel function"
import scipy.signal		# for sepfir2d() and gaussian()


class ScaleSpace_Blob :
	def __init__(self, idNum) :
		self.isBright = None
		self.signifVal = None
		self.appearance = None
		self.disappearance = None
		self.id_num = idNum

		# This will list grey blobs for scale each level.
		self.grey_blobs = []
		# This will list the support region for each scale level
		self.support_regions = []
		# This will list the scale value for each scale level
		self.scale_levels = []
		# This will list the bifurcations events for each scale level
		self.events = []

		self.approp_scalelvl = None
		self.approp_greyblob = None
		self.atBoundary = None


	def Start_ScaleBlob(self, greyBlob, scale_level) :
		self.appearance = scale_level
		self.grey_blobs.append(greyBlob)
		self.support_regions.append(ws.Support_Region(greyBlob.support_region))
		self.scale_levels.append(scale_level)

		newEvent = Scale_Event(Scale_Event.CREATION,
                                       [], [self],
                                       self.support_regions[-1].first_moment(),
                                       self.scale_levels[-1])
		self.events.append( [newEvent] )

		return newEvent

	def End_ScaleBlob(self) :
		self.disappearance = self.scale_levels[-1]
		endEvent = Scale_Event(Scale_Event.DESTRUCTION,
				       [self], [],
				       self.support_regions[-1].first_moment(),
				       self.scale_levels[-1])
		self.events.append( [endEvent] )

		return endEvent

	def Continue_ScaleBlob(self, greyBlob, scaleLevel) :
		self.support_regions.append(ws.Support_Region(greyBlob.support_region))
		self.grey_blobs.append(greyBlob)
		self.scale_levels.append(scaleLevel)

	def Split_ScaleBlob(self, greyBlob) :
		# So, we split a scale blob by updating the last event in the events list
		# Note that this function assumes that in order to know that a split has
		#      occurred, that the scale blob has already been "continued", and that
		#      the splits are being discovered one grey blob at a time.
		self.support_regions[-1].AddSupport(greyBlob.support_region)
		self.grey_blobs[-1].append(greyBlob)
		
		self.events[-1][-1].event_type = Scale_Event.SPLIT
		self.events[-1][-1].position = self.support_regions[-1].first_moment()

		return self.events[-1][-1]

#	def Merge_ScaleBlob(self, greyBlob, scaleBlobs, scaleLevel) :
#
#		self.events[-1][-1].event_type = Scale_Event.MERGE
#		self.events[-1][-1].scaleblobs_above += scaleBlobs
#		# TODO: Hmm, need to also update position, but how?
#
#		return self.events[-1][-1]
		

	def Add_Greylevel_Blobs(self, blobs, scale_lev) :
		new_support = ws.Support_Region()

		for aBlob in blobs :
			aBlob.scalespace_blob = self
			new_support.AddSupport(aBlob.support_region)

		self.support_regions.append(new_support)
		self.grey_blobs.append(blobs)
		self.scale_levels.append(scale_lev)

	def Add_Greylevel_Blob(self, aBlob) :
		# Assumes that we are adding one blob to the current level.
		# Therefore, a level must exist first...
		aBlob.scalespace_blob = self
		self.support_regions[-1].AddSupport(aBlob.support_region)
		self.grey_blobs[-1].append(aBlob)




class Scale_Event :
	CREATION = 0
	SPLIT = 1
	MERGE = 2
	DESTRUCTION = 3

	def __init__(self, event_type, scaleblobs_above = [], scaleblobs_below = [], position = None, scale_lev = None) :
		self.event_type = event_type
		self.scaleblobs_above = scaleblobs_above
		self.scaleblobs_below = scaleblobs_below
		self.position = position
		self.scale_level = scale_lev




class Scale_Level :
	def __init__(self, greyBlobs, image, greyMarks, scaleVal) :
		self.scaleVal = scaleVal
		self.image = image
		self.greyMarks = greyMarks
		self.greyBlobs = greyBlobs



class Primal_Sketch :
	UNMARKED = -1

	def __init__(self) :
		self.scale_levels = []

		self.scaleBlobs_bright = []
		self.events_bright = []

		self.scaleBlob_Marks = []




	def CreateSketch(self, image, scale_values) :
		prevScaleMarks = numpy.empty(image.shape, dtype=int)
		prevScaleMarks.fill(self.UNMARKED)

		for aScale in scale_values :
			if aScale == 0 :
				newImage = image.copy()
			else :
				newImage = self.DoConvolve(image, aScale, (4 * (aScale // 2)) + 3).astype(int)

			print "At level: ", aScale #, "  Image max:", newImage.max(), "   Image min:", newImage.min()
			
			greyblobs, greyMarks = ws.Watershed_Transform(newImage)

			newScale = Scale_Level(greyblobs, image, greyMarks, aScale)
			self.scale_levels.append(newScale)

			newScaleMarks = self.Link_Greyblobs(prevScaleMarks, newScale)
			self.scaleBlob_Marks.append(newScaleMarks)
			prevScaleMarks = newScaleMarks




	def DoConvolve(self, image, scale_level, winSize) :
		# NOTE: I know this isn't technically the best approach.
		# There is supposedly a better way using Bessel functions,
		#  but until I get a better idea how this is supposed to be
		#  implemented, I will do it this way.
		kernel = scipy.signal.gaussian(winSize, scale_level)
		return scipy.signal.sepfir2d(image, kernel / kernel.sum(),
						    kernel / kernel.sum())

#	def AddNewScaleLevel(self, image, greyblobs, greyMarks, scaleVal) :
#		# Right now, assume that the scale size is changing monotonically.
#		# Therefore, I won't bother with trying to sort and mess around with linkage issues.
#		self.scale_levels.append(Scale_Level(greyblobs, image, greyMarks, scaleVal))
		
	def Link_Greyblobs(self, prevScaleMarks, newScale) :
		currIDNum = len(self.scaleBlobs_bright)
		currScaleMarks = numpy.empty(prevScaleMarks.shape, dtype=int)
		currScaleMarks.fill(self.UNMARKED)

		ignoreThese = frozenset([self.UNMARKED])
		scaleBlobsMatched = set([])
		scaleMap = {}
		
		
		for aGreyBlob in newScale.greyBlobs :
			# Find out what scale space blob existed at the previous scale level at the location
			# of this grey blob's extremum. We automatically removed any UNMARKED as well.
			blobIndices = set([prevScaleMarks[anIndex] for anIndex in aGreyBlob.extremum]) - ignoreThese

			if len(blobIndices) == 0 :
				# This is an absolutely brand-new scale blob!
				for anIndex in aGreyBlob.support_region :
					currScaleMarks[anIndex] = currIDNum

				new_blob = ScaleSpace_Blob(currIDNum)
				new_event = new_blob.Start_ScaleBlob(aGreyBlob, newScale)

				currIDNum += 1
				self.scaleBlobs_bright.append(new_blob)
				self.events_bright.append(new_event)
			else :
				if len(blobIndices) > 1 :
					print "Degenerate situation! Not correctly implemented!"

				blobIndex = list(blobIndices)[0]	# just grabbing one
				# There is at least a continuation, but it could be many other things.
				if blobIndex in scaleBlobsMatched :
					scaleMap[blobIndex].append(aGreyBlob)
				else :
					scaleBlobsMatched.add(blobIndex)
					scaleMap[blobIndex] = [aGreyBlob]

		for (scaleBlobIndex, greyBlobs) in scaleMap.items() :
			aScaleBlob = self.scaleBlobs_bright[scaleBlobIndex]
			if len(greyBlobs) == 1 :
				# Is it a simple linkage or a merger?
				theSaddle = aScaleBlob.grey_blobs[-1].saddle
				if theSaddle is None or len(theSaddle.grey_blobs) <= 1 :
					# Simple linkage
					aScaleBlob.Continue_ScaleBlob(greyBlobs[0], newScale)
				else :
					# Merger
					currIDNum = self.Merge_ScaleBlobs(greyBlobs[0], currIDNum, currScaleMarks,
									  newScale, [aScaleBlob])

			else :
				# Is it a simple split or creation of new scale blobs?
				# NOTE: FOR NOW we are gonna lie!
				currIDNum = self.Split_ScaleBlob(greyBlobs, currIDNum, currScaleMarks,
								 newScale, aScaleBlob)
				



		# Any scale blobs from the previous level that are unmatched needs to be discontinued.
		unmatchedBlobs = set(prevScaleMarks.flat) - scaleBlobsMatched - ignoreThese
		for blobIndex in unmatchedBlobs :
			end_event = self.scaleBlobs_bright[blobIndex].End_ScaleBlob()
			self.events_bright.append(end_event)

		return currScaleMarks


	def Split_ScaleBlob(self, greyBlobs, currIDNum, currScaleMarks, scaleLevel, scaleBlob) :
		splitEvent = Scale_Event(Scale_Event.SPLIT,
                                         [scaleBlob], [],
                                         scaleBlob.support_regions[-1].first_moment(),
                                         scaleLevel)

		newScaleBlobs = []
		for aGreyBlob in greyBlobs :
			for anIndex in aGreyBlob.support_region :
				currScaleMarks[anIndex] = currIDNum

			newBlob = ScaleSpace_Blob(currIDNum)
			newBlob.grey_blobs.append(aGreyBlob)
			newBlob.appearance = scaleLevel
			newBlob.support_regions.append(aGreyBlob.support_region)
                	newBlob.scale_levels.append(scaleLevel)
			newBlob.events.append(splitEvent)

			currIDNum += 1

			newScaleBlobs.append(newBlob)

		# TODO: Probably some more things I was supposed to do...
		splitEvent.scaleBlobs_below = newScaleBlobs
		scaleBlob.events.append(splitEvent)
		self.events_bright.append(splitEvent)
		self.scaleBlobs_bright += newScaleBlobs

		return currIDNum



	def Merge_ScaleBlobs(self, greyBlob, currIDNum, currScaleMarks, scaleLevel, scaleBlobs) :
		for anIndex in greyBlob.support_region :
                	currScaleMarks[anIndex] = currIDNum
		
		newBlob = ScaleSpace_Blob(currIDNum)
		newBlob.grey_blobs.append(greyBlob)

		currIDNum += 1

		# Not exactly sure how we are going to represent it correctly,
		#     so we will get away with just creating a new blob for now
		
		newBlob.appearance = scaleLevel
		newBlob.support_regions.append(greyBlob.support_region)
                newBlob.scale_levels.append(scaleLevel)


		# I know, I am missing the other scale blobs...
		mergeEvent = Scale_Event(Scale_Event.MERGE,
					 scaleBlobs, [newBlob],
                                       	 greyBlob.support_region.first_moment(),
					 scaleLevel)

		newBlob.events.append(mergeEvent)
		for aScaleBlob in scaleBlobs :
			aScaleBlob.events.append(mergeEvent)
					
		# TODO: end the other scale blobs!
					

		self.events_bright.append(mergeEvent)
		self.scaleBlobs_bright.append(newBlob)

		return currIDNum
