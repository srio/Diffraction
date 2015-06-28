"""


fresnel_fourier_1D.py: calculates 1D fresnel diffraction via convolution by Fourier transform
 

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"

import numpy

if __name__ == '__main__':


    wavelength   =   5000e-10
    sourcesize   =   500e-6
    detsize = 8000e-6


    npoints = 100
    #detpoints =  1000
    distance =   1.00e0
    lensF        =   None
    
    position_x = numpy.linspace(-detsize/2,detsize/2,npoints)
    xoffset = position_x[0]
    xdelta = position_x[1] - position_x[0]
    xsize = npoints


    fields1 = numpy.ones(npoints) + 0j
    fields1[numpy.where(position_x < -0.5*sourcesize)] = 0.0
    fields1[numpy.where(position_x >  0.5*sourcesize)] = 0.0

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

    #back fft
    fields2 = numpy.fft.ifft(wfou_fft)

    fieldIntensity = numpy.abs(fields2)**2
    fieldPhase = numpy.arctan2(numpy.real(fields2), \
                               numpy.imag(fields2))
    #plots
    from matplotlib import pylab as plt

    plt.figure(1)
    plt.plot(position_x,fieldIntensity)

    plt.xlabel="Z [m]"
    plt.ylabel="Intensity [a.u.]"
    plt.title="Source "
    plt.show()

    #
    # write spec formatted file
    #
    out_file = "fresnel_fourier_1D.spec"
    f = open(out_file, 'w')
    header="#F %s \n\n#S  1 kirchhoff \n#N 3 \n#L Z[m]  intensity  phase\n"%out_file
    f.write(header)
    for i in range(len(position_x)):
        out = numpy.array((position_x[i],  fieldIntensity[i],fieldPhase[i]))
        f.write(("%20.11e "*out.size+"\n") % tuple( out.tolist()))
    f.close()
    print ("File written to disk: %s"%out_file)
