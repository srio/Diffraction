#
# plots radial dependency of Airy disk
#
# (c) srio@esrf.eu 20150624
#

from scipy.special import jv
import numpy

#
#inputs (in SI units)
#
# wavelength = 1e-10 # 500e-9 # 1e-10
# distance = 1.0
# aperture_diameter = 1e-6 # 1e-3 # 1e-6

wavelength = 5000e-10 # 500e-9 # 1e-10
distance = 1.0
aperture_diameter = 500e-6 # 1e-3 # 1e-6


sin_theta = 1.22*wavelength/aperture_diameter
print("Angular radius of first Airy ring: %15.10f"%sin_theta)
print("spatial position at first minimum: %15.10f mm"%(sin_theta*distance*1e3))

sin_theta_array = numpy.linspace(-3*sin_theta,3*sin_theta,1000)
x = (2*numpy.pi/wavelength) * (aperture_diameter/2) * sin_theta_array 
x_over_pi = x / numpy.pi
U_vs_theta = 2*jv(1,x)/x
I_vs_theta = U_vs_theta**2

#CALCULATE fwhm
tt = numpy.where(I_vs_theta>=max(I_vs_theta)*0.5)
if I_vs_theta[tt].size > 1:
    binSize = sin_theta_array[1]-sin_theta_array[0]
    FWHM = binSize*(tt[0][-1]-tt[0][0])
print("lambda/(Dx Dtheta): %15.10f"%(wavelength/(aperture_diameter*2*sin_theta)))
print("lambda/(FWHMx FWHMtheta): %15.10f"%(wavelength/(aperture_diameter*FWHM)))
print("lambda/(SIGMAx SIGMAtheta): %15.10f"%(wavelength/(aperture_diameter*FWHM/2.35/2.35)))

#
# range of validity
#
print("Fraunhoffer diffraction valid for distances > > a^2/lambda = %f m"%((aperture_diameter/2)**2/wavelength))
#
#
#plots
#
from matplotlib import pylab as plt

plt.figure(1)
plt.plot(x_over_pi,I_vs_theta)
plt.xlabel("k a sin(theta) / pi")
plt.ylabel("I/Io")

plt.figure(2)
plt.plot(sin_theta_array,I_vs_theta)
plt.xlabel("theta")
plt.ylabel("I/Io")

plt.figure(3)
plt.plot(sin_theta_array*distance,I_vs_theta)
plt.xlabel("x [m]")
plt.ylabel("I/Io")
plt.show()


