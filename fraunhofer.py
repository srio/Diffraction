"""


fraunhofer.py: calculates 2D Fraunhofer diffraction  (via Fourier Transform)


"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"


import numpy as np
 
#
# inputs (in SI)
#

wavelength        = 1.24e-10
aperture_diameter = 40e-6 # if Gaussian, aperture_diameter = 2.35*sigma

pixelsize         = 1e-6
npixels           = 1024

aperture_type     = 0    # 0=circular, 1=Gaussian


#
#create array at object (aperture) plane
#
p_i = np.arange(1,npixels)
p_x = p_i - npixels/2
p_y = npixels/2 - p_i


p_xx = p_x[:, np.newaxis]
p_yy = p_y[np.newaxis, :]

image = np.zeros((npixels-1,npixels-1))  # amplitude map

if aperture_type == 0:  # Circular aperture
    radius = (aperture_diameter/2)/pixelsize # in pixels
    print("radius=%f pixels"%radius)
    image_illuminated_indices = np.where(p_xx**2 + p_yy**2 < radius**2)
    image[image_illuminated_indices] = 1.0
elif aperture_type == 1:  # Gaussian
    sigma = aperture_diameter/pixelsize/2.35
    print("sigma=%f pixels"%sigma)
    rho2 = p_xx**2 + p_yy**2
    #TODO: add Gaussian amplitude
    image = np.exp(-rho2/2/sigma**2) # amplitude
else:
    raise ValueError("Aperture type (shape) not valid")


#
#compute Fourier transform 
#

F1 = np.fft.fft2(image)  # Take the fourier transform of the image.
# Now shift the quadrants around so that low spatial frequencies are in
# the center of the 2D fourier transformed image.
F2 = np.fft.fftshift( F1 )

# abscissas
freq_nyquist = 0.5/pixelsize
freq_n = np.linspace(-1.0,1.0,npixels-1) 
freq = freq_n * freq_nyquist
freq = freq * wavelength

ka = 2*np.pi/wavelength * (aperture_diameter/2)
x = freq*ka
x_over_pi = x / np.pi
 

#
# range of validity
#
print("Fraunhoffer diffraction valid for distances > > a^2/lambda = %f m"%((aperture_diameter/2)**2/wavelength))

 

#
#CALCULATE fwhm 's
#
IntensityH = np.abs(F2[int(npixels/2),:])**2 # horizontal profile
tt = np.where(IntensityH>=max(IntensityH)*0.5)
if IntensityH[tt].size > 1:
    binSize = freq[1]-freq[0]
    FWHM_theta = binSize*(tt[0][-1]-tt[0][0])
    print("FWHM of intensity pattern: %f rad"%FWHM_theta)

    ImageH = np.abs(image[int(npixels/2),:])**2
    tt = np.where(ImageH>=max(ImageH)*0.5)
    if ImageH[tt].size > 1:
        FWHM_x = pixelsize*(tt[0][-1]-tt[0][0])
        print("FWHM of original aperture: %f m"%FWHM_x)
        print("lambda/(FWHMx FWHMtheta): %15.10f"%(wavelength/(FWHM_theta*FWHM_x)))
        print("lambda/(SIGMAx SIGMAtheta): %15.10f"%(wavelength/(FWHM_theta*FWHM_x/2.35/2.35)))



#
#make plots
#

import matplotlib.pylab as plt

# Now plot up both
plt.figure(1)
plt.clf()
plt.imshow(  image , extent=[1e6*freq[0],1e6*freq[-1],1e6*freq[0],1e6*freq[-1]], cmap=plt.cm.Greys)
plt.title("aperture intensity")
plt.xlabel("X [um]")
plt.ylabel("Z [um]")
plt.colorbar()

plt.figure(2)
plt.clf()
plt.imshow( np.abs(F2), extent=[1e6*freq[0],1e6*freq[-1],1e6*freq[0],1e6*freq[-1]], )
plt.title("diffracted intensity")
plt.xlabel("X [um]")
plt.ylabel("Z [um]")
plt.colorbar()

plt.figure(3)
plt.clf()
if aperture_type == 0: #circular, also display analytical values
    from scipy.special import jv
    x = (2*np.pi/wavelength) * (aperture_diameter/2) * freq
    U_vs_theta = 2*jv(1,x)/x 
    I_vs_theta = U_vs_theta**2 * IntensityH.max()
    plt.plot( freq*1e6, IntensityH, '-', freq*1e6 , I_vs_theta, '.' )
elif aperture_type == 1: #Gaussian
    plt.plot( freq*1e6, IntensityH)

plt.title("Horizontal profile of diffracted intensity")
plt.xlabel('theta [urad]')
plt.ylabel('Diffracted intensity [a.u.]')
 
plt.show()
