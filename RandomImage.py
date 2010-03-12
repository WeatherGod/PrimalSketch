import numpy


def RandomImage(pointCnt, maxRadius, maxPower, domainShape) :
	xPos = numpy.random.randint(0, domainShape[0], pointCnt)
	yPos = numpy.random.randint(0, domainShape[1], pointCnt)
	theRadius = numpy.random.randint(1, maxRadius + 1, pointCnt)
	thePower = numpy.random.rand(pointCnt) * maxPower

	dataVals = numpy.zeros(domainShape)

	for x in range(domainShape[0]) :
		for y in range(domainShape[1]) :
			dists = numpy.sqrt((x - xPos)**2 + (y - yPos)**2)
			dataVals[x, y] = numpy.sum(thePower * numpy.exp(-(dists**2) / (2*(theRadius/2.0)**2)))

	return dataVals
			
