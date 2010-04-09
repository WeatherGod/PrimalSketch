#!/usr/bin/env python

from LoadRastRadar import *
import ScaleSpace as ss
import numpy

import pylab
import mpl_toolkits.mplot3d.axes3d as p3

scales = [1, 5, 7, 11, 21, 31, 41]

def foo(ps) :
	radarData = LoadRastRadar("6500KTLX20050514_052255.nc")
	radarData['vals'] = numpy.nan_to_num(numpy.squeeze(radarData['vals']))

	print "Data min:", radarData['vals'].min(), "   Data max:", radarData['vals'].max()

	imageData = (radarData['vals'] - radarData['vals'].min()).astype(int)

	print "Image max: ", imageData.max()


	ps.CreateSketch(imageData, scales)

def Make_ScaleSpace_Plot(ps) :
	someColors = ['r', 'g', 'm', 'b', 'c', 'k']
	X, Y = numpy.meshgrid(numpy.arange(1200), numpy.arange(925))


	fig = pylab.figure()
	ax = p3.Axes3D(fig)
	pylab.hold(True)

	for index, aBlob in enumerate(ps.scaleBlobs_bright) :
		
		Z = numpy.zeros(Y.shape)
		(y, x, z) = Blob2Coords(aBlob)
		Z[y, x] = z
		ax.contour3D(X, Y, Z, colors = someColors[index % 6])


	pylab.show()


def Blob2Coords(blob) :
        x = []
        y = []
        z = []

        for aRegion, scaleVal in zip(blob.support_regions, blob.scale_values) :
                xtemp, ytemp = zip(*list(aRegion.pixels))
                x += xtemp
                y += ytemp
                z += [scaleVal + 1] * len(xtemp)

        return x, y, z
		

