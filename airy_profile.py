"""


airy_profile.py: calculates 1D radial intensity of the Fraunhofer diffraction of a disk aperture  (Airy pattern)


"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"


from scipy.special import jv
import numpy

#
#inputs (in SI units)
#

# wavelength = 1e-10 # 500e-9 # 1e-10
# distance = 1.0
# aperture_diameter = 1e-6 # 1e-3 # 1e-6

# wavelength = 5000e-10 # 500e-9 # 1e-10
# distance = 1.0
# aperture_diameter = 500e-6 # 1e-3 # 1e-6

wavelength = 1.24e-10 # 10keV
distance = 3.6
aperture_diameter = 40e-6 # 1e-3 # 1e-6


# sin(theta) ~ theta

sin_theta = 1.22*wavelength/aperture_diameter

print("Angular radius of first Airy ring: %15.10f rad"%sin_theta)
print("spatial position at first minimum: %15.10f mm"%(sin_theta*distance*1e3))

sin_theta_array = numpy.linspace(-3*sin_theta,3*sin_theta,1000)
x = (2*numpy.pi/wavelength) * (aperture_diameter/2) * sin_theta_array 
x_over_pi = x / numpy.pi
electric_field = 2*jv(1,x)/x
intensity = electric_field**2

#
#CALCULATE fwhm
#
tt = numpy.where(intensity>=max(intensity)*0.5)
if intensity[tt].size > 1:
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
# write spec formatted file
#
out_file = "airy_profile.spec"
f = open(out_file, 'w')
header="#F %s \n\n#S  1 airy profile \n#N 4 \n#L ka sin(theta)/pi  sin(theta)  Z[m]  intensity\n"%out_file
f.write(header)
for i in range(len(x_over_pi)):
    out = numpy.array((x_over_pi[i],  sin_theta_array[i],sin_theta_array[i]*distance, intensity[i]))
    f.write(("%20.11e "*out.size+"\n") % tuple( out.tolist()))
f.close()
print ("File written to disk: %s"%out_file)

#
#
#plots
#
from matplotlib import pylab as plt

plt.figure(1)
plt.plot(x_over_pi,intensity)
plt.xlabel("k a sin(theta) / pi")
plt.ylabel("I/Io")
plt.title("Airy pattern profile")

plt.figure(2)
plt.plot(sin_theta_array,intensity)
plt.xlabel("theta ~ sin(theta)")
plt.ylabel("I/Io")
plt.title("Airy pattern profile")

plt.figure(3)
plt.plot(sin_theta_array*distance,intensity)
plt.xlabel("x [m]")
plt.ylabel("I/Io")
plt.title("Airy pattern profile")

plt.show()