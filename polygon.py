from PyMca5.PyMca import ctools
pnpoly = ctools.pnpoly
import numpy

vertices = numpy.zeros(8)
vertices.shape = 4, 2
vertices[0,:] = [0.0, 0.0]
vertices[1,:] = [1.0, 0.0]
vertices[2,:] = [1.0, 1.0]
vertices[3,:] = [0.0, 1.0]

print(" 0.5, 0.5 is inside ", pnpoly(vertices, [[0.5, 0.5]]))
print(" 1.5, 0.5 is inside ", pnpoly(vertices, [[1.5, 0.5]]))
