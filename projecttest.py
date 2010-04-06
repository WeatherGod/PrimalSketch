#!/usr/bin/env python

from LoadRastRadar import *
import ScaleSpace as ss
import numpy

def foo() :
	radarData = LoadRastRadar("6500KTLX20050514_052255.nc")
	radarData['vals'] = numpy.nan_to_num(numpy.squeeze(radarData['vals']))

	print "Data min:", radarData['vals'].min(), "   Data max:", radarData['vals'].max()

	imageData = (radarData['vals'] - radarData['vals'].min()).astype(int)

	print "Image max: ", imageData.max()


	ps = ss.Primal_Sketch()

	ps.CreateSketch(imageData, [])


