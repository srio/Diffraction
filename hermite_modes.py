"""


hermite_modes.py: calculates Gauss-Schell modes

see e.g. https://commons.wikimedia.org/wiki/File:Hermite-gaussian.png

"""

__author__ = "M. Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2016"



import numpy
from scipy.special import hermite
from srxraylib.plot.gol import plot_image


#
# inputs
#
size_x = 1.0
size_y = .50
n_x = 620
n_y = 320
w = size_x/8.
m = 3
n = 0

X = numpy.linspace(-size_x/2,size_x/2,n_x)
Y = numpy.linspace(-size_y/2,size_y/2,n_y)

XX = numpy.outer(X,numpy.ones_like(Y))
YY = numpy.outer(numpy.ones_like(X),Y)

out =     (hermite(m)(numpy.sqrt(2)*XX/w)*numpy.exp(-XX**2/w**2))**2 \
        * (hermite(n)(numpy.sqrt(2)*YY/w)*numpy.exp(-YY**2/w**2))**2

plot_image(out,X,Y)
