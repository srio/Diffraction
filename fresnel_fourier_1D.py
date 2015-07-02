"""


fresnel_fourier_1D.py: calculates 1D fresnel diffraction via convolution by Fourier transform


"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"


import numpy

if __name__ == '__main__':


    # wavelength   =   5000e-10
    # aperture_diameter   =   500e-6
    # detector_size = 8000e-6
    # distance =   1.00e0



    wavelength = 1.24e-10 # 10keV
    aperture_diameter = 40e-6 # 1e-3 # 1e-6
    detector_size = 800e-6
    distance = 3.6

    npoints = 1000

    lensF        =   None
    
    position_x = numpy.linspace(-detector_size/2,detector_size/2,npoints)
    xoffset = position_x[0]
    xdelta = position_x[1] - position_x[0]
    xsize = npoints


    fields1 = numpy.ones(npoints) + 0j
    fields1[numpy.where(position_x < -0.5*aperture_diameter)] = 0.0
    fields1[numpy.where(position_x >  0.5*aperture_diameter)] = 0.0

    #apply ideal lens
    if 0:
        focallength = 100.0
        knum = 2.0*numpy.pi/wavelength
        fields1 *= numpy.exp(-knum*1j*position_x**2/focallength/2.0)

    fft_size = npoints
    fft_delta = 1.0/xsize/xdelta
    if numpy.mod(npoints,2) == 1:
        fft_offset = -fft_delta*float(npoints-1)/2.0
    else:
        fft_offset = -fft_delta*float(npoints)/2.0

    # FT
    F1 = numpy.fft.fft(fields1)
    wfou_fft = numpy.fft.fftshift(F1)
    wfou_fft_x = numpy.arange(start=fft_offset, stop = -fft_offset, step=fft_delta, )

    #propagate
    wfou_fft *= numpy.exp(-1j * numpy.pi * wavelength * distance * wfou_fft_x**2 )

    #back FT
    fields2 = numpy.fft.ifft(wfou_fft)

    fieldIntensity = numpy.abs(fields2)**2
    fieldPhase = numpy.arctan2(numpy.real(fields2), \
                               numpy.imag(fields2))

    #
    # write spec formatted file
    #
    out_file = "fresnel_fourier_1D.spec"
    f = open(out_file, 'w')
    header="#F %s \n\n#S  1 fresnel diffraction \n#N 3 \n#L X[m]  intensity  phase\n"%out_file
    f.write(header)
    for i in range(len(position_x)):
        out = numpy.array((position_x[i],  fieldIntensity[i],fieldPhase[i]))
        f.write(("%20.11e "*out.size+"\n") % tuple( out.tolist()))
    f.close()
    print ("File written to disk: %s"%out_file)

    #plots


    from matplotlib import pylab as plt
    from scipy.special import jv
    sin_theta = position_x / distance
    x = (2*numpy.pi/wavelength) * (aperture_diameter/2) * sin_theta
    U_vs_theta = 2*jv(1,x)/x
    I_vs_theta = U_vs_theta**2 * fieldIntensity.max()


    plt.figure(1)

    plt.plot(position_x*1e6,fieldIntensity,'-',position_x*1e6,I_vs_theta,)
    plt.title("Fresnel Diffraction")
    plt.xlabel("X [um]")
    plt.ylabel("Intensity [a.u.]")
    plt.show()


