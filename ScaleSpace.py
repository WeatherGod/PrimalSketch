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
		self.grey_blobs = []
		self.approp_scalelvl = None
		self.approp_greyblob = None
		self.atBoundary = None


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
		print "Base scale"
		greyblobs, basinMarks = ws.Watershed_Transform(image)
		AddBlobs(greyblobs, basinMarks, 0.0)

		for aScale in scale_levels :
			print "At level: ", aScale
			greyblobs, basinMarks = ws.Watershed_Transform(DoConvolve(image, aScale, (4 * (aScale // 2)) + 1))
			AddBlobs(greyblobs, basinMarks, aScale)




	def DoConvolve(self, image, scale_level, winSize) :
		# NOTE: I know this isn't technically the best approach.
		# There is supposedly a better way using Bessel functions,
		#  but until I get a better idea how this is supposed to be
		#  implemented, I will do it this way.
		return scipy.signal.sepfir2d(image, scipy.signal.gaussian(winSize, scale_level),
						    scipy.signal.gaussian(winSize, scale_level))

	def AddBlobs(self, greyblobs, image, scalesize) :
		self.scale_levels.append(Scale_Level(greyblobs, image, scalesize))


