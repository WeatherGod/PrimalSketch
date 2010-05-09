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
		self.grey_blobs.append([greyBlob])
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
		self.grey_blobs.append([greyBlob])
		self.scale_levels.append(scaleLevel)

		scaleEvent = Scale_Event(Scale_Event.CONTINUE,
                                         [self], [self],
                                         self.support_regions[-1].first_moment(),
                                         self.scale_levels[-1])
		
		self.events.append( [scaleEvent] )

		return scaleEvent

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

	def Merge_ScaleBlob(self, greyBlob, scaleBlobs, scaleLevel = None) :
		# Well, it doesn't really "merge" yet, but it does update some of the info.
		if scaleLevel is not None :
			self.Continue_ScaleBlob(greyBlob, scaleLevel)

		self.events[-1][-1].event_type = Scale_Event.MERGE
		self.events[-1][-1].scaleblobs_above += scaleBlobs
		# TODO: Hmm, need to also update position, but how?

		return self.events[-1][-1]
		

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
	CONTINUE = 2
	MERGE = 3
	DESTRUCTION = 4

	def __init__(self, event_type, scaleblobs_above = [], scaleblobs_below = [], position = None, scale_lev = None) :
		self.event_type = event_type
		self.scaleblobs_above = scaleblobs_above
		self.scaleblobs_below = scaleblobs_below
		self.position = position
		self.scale_level = scale_lev




class Scale_Level :
	def __init__(self, greyBlobs, image, blobMarks, scaleVal) :
		self.scaleVal = scaleVal
		#self.image = image
		self.blobMarks = blobMarks
		self.greyBlobs = greyBlobs



class Primal_Sketch :
	UNMARKED = -1

	def __init__(self) :
		self.scale_levels = []

		self.scaleBlobs_bright = []
		self.events_bright = []

		self.scaleBlob_Marks = []




	def CreateSketch(self, image, scale_values) :

		for aScale in scale_values :
			if aScale == 0 :
				# Maybe this should be a copy?
				newImage = image
			else :
				newImage = self.DoConvolve(image, aScale, (4 * (aScale // 2)) + 3).astype(int)

			print "At level: ", aScale #, "  Image max:", newImage.max(), "   Image min:", newImage.min()
			
			greyblobs, blobMarks = ws.Watershed_Transform(newImage)
			self.AddNewScaleLevel(newImage, greyblobs, blobMarks, aScale)
			self.Create_ScaleSpace_Blobs()




	def DoConvolve(self, image, scale_level, winSize) :
		# NOTE: I know this isn't technically the best approach.
		# There is supposedly a better way using Bessel functions,
		#  but until I get a better idea how this is supposed to be
		#  implemented, I will do it this way.
		kernel = scipy.signal.gaussian(winSize, scale_level)
		return scipy.signal.sepfir2d(image, kernel / kernel.sum(),
						    kernel / kernel.sum())

	def AddNewScaleLevel(self, image, greyblobs, greyMarks, scaleVal) :
		# Right now, assume that the scale size is changing monotonically.
		# Therefore, I won't bother with trying to sort and mess around with linkage issues.
		self.scale_levels.append(Scale_Level(greyblobs, image, greyMarks, scaleVal))
		newMarks = numpy.empty(greyMarks.shape, dtype=int)
		newMarks.fill(self.UNMARKED)
		self.scaleBlob_Marks.append(newMarks)
		

	def Create_ScaleSpace_Blobs(self) :
		# Skip out early if there is nothing to process
		if len(self.scale_levels) == 0 : return

		currIDNum = len(self.scaleBlobs_bright)

		atLevel = self.scale_levels[-1]
		
		if len(self.scale_levels) > 1 :
			previousBlobMarks = self.scaleBlob_Marks[-2]
		else :
			# if there is no 'previous level', then make a fake one.
			previousBlobMarks = numpy.empty(atLevel.blobMarks.shape, dtype=int)
			previousBlobMarks.fill(self.UNMARKED)
		
		ignoreThese = frozenset([self.UNMARKED])
		scaleBlobsMatched = set([])

		for aGreyBlob in atLevel.greyBlobs :
			blobIndices = set([previousBlobMarks[anIndex] for anIndex in aGreyBlob.support_region]) - ignoreThese

			# So, if len(blobIndices) is zero, then we have a new scalespace blob!
			#     if len(blobIndices) is one, then we have a continuation and possibly a split
			#     if len(blobIndices) is greater than one, we have a merge, and maybe a split as well

			if len(blobIndices) == 0 :

				for anIndex in aGreyBlob.support_region :
					self.scaleBlob_Marks[-1][anIndex] = currIDNum

				new_blob = ScaleSpace_Blob(currIDNum)
				new_event = new_blob.Start_ScaleBlob(aGreyBlob, atLevel)

				currIDNum += 1
				self.scaleBlobs_bright.append(new_blob)
				self.events_bright.append(new_event)

			elif len(blobIndices) == 1 :
				blobIndex = list(blobIndices)[0]
				print blobIndex
				for anIndex in aGreyBlob.support_region :
					self.scaleBlob_Marks[-1][anIndex] = blobIndex

				# If a scale blob is matched at least twice, then we have a split,
				if blobIndex in scaleBlobsMatched :
					self.scaleBlobs_bright[blobIndex].Split_ScaleBlob(aGreyBlob)
				else :
					contEvent = self.scaleBlobs_bright[blobIndex].Continue_ScaleBlob(aGreyBlob, atLevel)
					self.events_bright.append(contEvent)
					scaleBlobsMatched.add(blobIndex)

			else :
				# TODO: Figure out what to do here...
				print "Grrr... a bad component labeling issue.  Just pick the first one and go with it."
				# Might be a smarter way to choose...  maybe choose one that already has been matched?
				blobIndex = list(blobIndices)[0]

				for anIndex in aGreyBlob.support_region :
					self.scaleBlob_Marks[-1][anIndex] = blobIndex

				# If a scale blob is matched at least twice, then we have a split,
				# but is this possible when we already have a merger?
				# I don't think it is possible and there might be some geometric proof saying that
				# such an event isn't a real split/merge.  Looks like a paper is saying that the
				# event is actually more primitive events at unresolved scale levels.
                                if blobIndex in scaleBlobsMatched :
					print "ARRGH!  This isn't supposed to happen!  Don't worry, it isn't a critical issue right now."
					# TODO: We are just gonna have to get away with incorrect information for now...
                                        self.scaleBlobs_bright[blobIndex].Add_Greylevel_Blob(aGreyBlob)
                                else :
					mergeEvent = self.scaleBlobs_bright[blobIndex].Merge_ScaleBlob(aGreyBlob, 
											       [self.scaleBlobs_bright[index]
												 for index in blobIndices],
												atLevel)
                                        #self.scaleBlobs_bright[blobIndex].Add_Greylevel_Blobs([aGreyBlob], atLevel)
					self.events_bright.append(mergeEvent)
                                        scaleBlobsMatched.add(blobIndex)
					


		# Any scale blobs from the previous level that are unmatched needs to be discontinued.
		unmatchedBlobs = set(previousBlobMarks.flat) - scaleBlobsMatched - ignoreThese
		for blobIndex in unmatchedBlobs :
			end_event = self.scaleBlobs_bright[blobIndex].End_ScaleBlob()
			self.events_bright.append(end_event)
			

