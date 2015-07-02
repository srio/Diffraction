"""

fresnel: 

        functions: 
             goFromTo: calculates the phase shift matrix
 

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2012"

import numpy, math

def goFromTo(source,image,distance=1.0,lensF=None,wavelength=1e-10):
    distance = numpy.array(distance)
    x1 = numpy.outer(source,numpy.ones(image.size))
    x2 = numpy.outer(numpy.ones(source.size),image)
    r = numpy.sqrt( numpy.power(x1-x2,2) + numpy.power(distance,2) )
    # add lens at the image plane
    if lensF != None:
      r = r - numpy.power(x1-x2,2)/lensF
    wavenumber = numpy.pi*2/wavelength
    return numpy.exp(1.j * wavenumber *  r)

if __name__ == '__main__':

    # wavelength   =   1e-10
    # aperture_diameter   =   10e-6
    # detector_size = 0.8e-3
    # #wavelength   =   500e-9
    # #aperture_diameter   =   1e-3
    # #detector_size = 4e-3
    #
    # sourcepoints = 1000
    # detpoints =  1000
    # distance =   1.00
    # lensF        =   None

    # wavelength   =   5000e-10
    # sourcesize   =   500e-6
    # detector_size = 0.008
    #wavelength   =   500e-9
    #aperture_diameter   =   1e-3
    #detector_size = 4e-3

    wavelength = 1.24e-10 # 10keV
    aperture_diameter = 40e-6 # 1e-3 # 1e-6
    detector_size = 800e-6
    distance = 3.6


    sourcepoints = 1000
    detpoints =  1000
    lensF        =   None

    sourcesize = aperture_diameter
    
    position1x = numpy.linspace(-sourcesize/2,sourcesize/2,sourcepoints)
    position2x = numpy.linspace(-detector_size/2,detector_size/2,detpoints)
    
    fields12 = goFromTo(position1x,position2x,distance, \
        lensF=lensF,wavelength=wavelength)
    print ("Shape of fields12: ",fields12.shape)

    #prepare results
    fieldComplexAmplitude = numpy.dot(numpy.ones(sourcepoints),fields12)
    print ("Shape of Complex U: ",fieldComplexAmplitude.shape)
    print ("Shape of position1x: ",position1x.shape)
    fieldIntensity = numpy.power(numpy.abs(fieldComplexAmplitude),2)
    fieldPhase = numpy.arctan2(numpy.real(fieldComplexAmplitude), \
                               numpy.imag(fieldComplexAmplitude))


    #
    # write spec formatted file
    #
    out_file = "fresnel_kirchhoff_1D.spec"
    f = open(out_file, 'w')
    header="#F %s \n\n#S  1 fresnel-kirchhoff diffraction integral\n#N 3 \n#L X[m]  intensity  phase\n"%out_file

    f.write(header)
    
    for i in range(detpoints):
       out = numpy.array((position2x[i], fieldIntensity[i], fieldPhase[i]))
       f.write( ("%20.11e "*out.size+"\n") % tuple( out.tolist())  )
    
    f.close()
    print ("File written to disk: %s"%out_file)

    #
    #plots
    #
    from matplotlib import pylab as plt

    plt.figure(1)
    plt.plot(position2x*1e6,fieldIntensity)
    plt.title("Fresnel-Kirchhoff Diffraction")
    plt.xlabel("X [um]")
    plt.ylabel("Intensity [a.u.]")
    plt.show()