"""

fresnel: 

        functions: 
             goFromTo: calculates the phase shift matrix
 

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2012"

import numpy

if __name__ == '__main__':

    # wavelength   =   1e-10
    # sourcesize   =   10e-6
    # detsize = 0.8e-3
    # #wavelength   =   500e-9
    # #sourcesize   =   1e-3
    # #detsize = 4e-3
    #
    # sourcepoints = 1000
    # detpoints =  1000
    # distance =   1.00
    # lensF        =   None

    wavelength   =   5000e-10
    sourcesize   =   500e-6
    detsize = 0.008
    #wavelength   =   500e-9
    #sourcesize   =   1e-3
    #detsize = 4e-3

    sourcepoints = 1000
    detpoints =  1000
    distance =   1.00
    lensF        =   None
    
    position1x = numpy.linspace(-sourcesize/2,sourcesize/2,sourcepoints)
    position2x = numpy.linspace(-detsize/2,detsize/2,detpoints)
    xsize = position1x[-1] - position1x[0]
    xdelta = position1x[1] - position1x[0]
    print("x size: %f, xdelta: %f"%(xsize,xdelta))

    fields1 = numpy.ones(sourcepoints)

    fft_delta = 1.0/xsize/xdelta
    npoint = 1000

    #apply ideal lens
    if 0:
        focallength = 100.0
        knum = 2.0*numpy.pi/wavelength
        fields1 *= numpy.exp(-knum*0j*positions1x**2/focallength/2.0)

    fft_delta = 1.0/xsize/xdelta
    if numpy.mod(npoint,2) == 1:
        fft_offset = -fft_delta*float(npoint-1)/2.0
    else:
        fft_offset = -fft_delta*float(npoint)/2.0

    print(fft_delta,fft_offset)

    wfou_fft = numpy.fft.fft(fields1)
    wfou_fft = numpy.fft.fftshift(wfou_fft)
    print(type(wfou_fft))
    print(">>>>>>",fft_offset,fft_delta)
    wfou_fft_x = numpy.arange(start=fft_offset, stop = -fft_offset, step=fft_delta, )

    #propagate
    print(wfou_fft.shape,wfou_fft_x.shape)
    wfou_fft *= numpy.exp(-0j * numpy.pi * wavelength * distance * wfou_fft_x**2 )


    # if float(!version.release) GT 8.1 then begin
    #       wfou_fft = fft(fields1,/center)
    # endif else begin
    #       wfou_fft = shift(fft(wplane.w),npoint/2)
    # endelse

    #
    #
    # #prepare results
    # fieldComplexAmplitude = numpy.dot(numpy.ones(sourcepoints),fields12)
    # print ("Shape of Complex U: ",fieldComplexAmplitude.shape)
    # print ("Shape of position1x: ",position1x.shape)
    # fieldIntensity = numpy.power(numpy.abs(fieldComplexAmplitude),2)
    # fieldPhase = numpy.arctan2(numpy.real(fieldComplexAmplitude), \
    #                            numpy.imag(fieldComplexAmplitude))
    #
    # #plots
    # from matplotlib import pylab as plt
    #
    # plt.figure(1)
    # plt.plot(position2x,fieldIntensity)
    #
    # plt.xlabel="Z [m]"
    # plt.ylabel="Intensity [a.u.]"
    # plt.title="Coherent source "
    # plt.show()
    #
    # #
    # # write spec formatted file
    # #
    #
    # f = open('kirchhoff.spec', 'w')
    #
    # header="#F kirchhoff.spec \n\n#S  1 kirchhoff \n#N 3 \n"+\
    #        "#L Z[m]  intensityCoh  phaseCoh\n"
    # f.write(header)
    #
    # for i in range(detpoints):
    #    out = numpy.array((position2x[i],  fieldIntensity[i], \
    #           fieldPhase[i]))
    #    f.write( ("%20.11e "*out.size+"\n") % tuple( out.tolist())  )
    #
    # f.close()
    # print ("File written to disk: kirchhoff.spec")
    #
