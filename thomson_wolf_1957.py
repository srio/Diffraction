
"""
Program: 
ThomsonWolf1957: Compute Intensity profiles in Fig. 4, Thompson and Wolf
                JOSA 47 (1957) 895-902

Functions: 
                goFromTo: calculates the phase shift matrix

"""
__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2012"

import numpy

def goFromTo(source,image,distance=1.0,lensF=None,wavelength=1e-10):
    #distance = numpy.array(distance)
    x1 = numpy.outer(source,numpy.ones(image.size))
    x2 = numpy.outer(numpy.ones(source.size),image)
    r = numpy.sqrt( numpy.power(x1-x2,2) + numpy.power(distance,2) )
    # add lens at the image plane
    if lensF != None:
      x10 = numpy.outer(source*0,numpy.ones(image.size))
      r = r - numpy.power(x10-x2,2)/lensF
    wavenumber = numpy.pi*2/wavelength
    return numpy.exp(1.j * wavenumber *  r)

if __name__ == '__main__':
    # inputs (in SI units)
    sourcesize   = 90e-6
    sourcepoints = 200
    slitdistance = 1.59
    slitgap      = numpy.array((6,8,10,12,17,23,25))*1e-3  # 2h
    slitsize     = 1400e-6  #
    slitpoints   = 500
    detdistance  = 1.59
    detsize      = 6e-3
    detpoints    = 3000
    wavelength   = 579.0e-9
    realisations = 1000
    lensF        = detdistance

    # open output file
    outfile='thomson_wolf_1957.spec'
    f = open(outfile, 'w')
    header="#F %s \n"%outfile
    f.write(header)

    #loop over slitgap
    for j in range(slitgap.size):
        #compute position arrays
        print ("Calculation for 2h=",slitgap[j])
        position1x = numpy.linspace(-sourcesize/2,sourcesize/2,sourcepoints)
        tmp = numpy.linspace(-slitsize/2,slitsize/2,slitpoints)
        if slitgap[j] != 0:
           position2x = numpy.concatenate((tmp-slitgap[j]/2,tmp+slitgap[j]/2))
        else:
           position2x = tmp
        position3x = numpy.linspace(-detsize/2,detsize/2,detpoints)
        
        fields12 = goFromTo(position1x,position2x,slitdistance, lensF=lensF,wavelength=wavelength)
        fields23 = goFromTo(position2x,position3x,detdistance, lensF=None,wavelength=wavelength)
        # from 1 to 3, matrix multiplication
        fields13 = numpy.dot(fields12,fields23)
    
        fieldComplexAmplitude = numpy.dot(numpy.ones(sourcepoints),fields13)
        fieldIntensity = numpy.power(numpy.abs(fieldComplexAmplitude),2)
        fieldPhase = numpy.arctan2(numpy.real(fieldComplexAmplitude), numpy.imag(fieldComplexAmplitude))
        
        # do the ensemble average
        tmpSource = numpy.exp(1.j*2*numpy.pi* numpy.random.mtrand.rand(sourcepoints))
        fieldSource=tmpSource
        fieldIntensityEA = numpy.power(numpy.abs(fieldComplexAmplitude),2)
        for i in range(realisations-1): 
          tmpSource = numpy.exp(1.j*2* numpy.pi*numpy.random.mtrand.rand(sourcepoints))
          fieldComplexAmplitude = numpy.dot( tmpSource, fields13)
          fieldIntensityEA = fieldIntensityEA + numpy.power(numpy.abs(fieldComplexAmplitude),2)
        
        header="\n#S "+str(j+1)+" 2h="+str(slitgap[j]*1e3)+" mm\n"
        f.write(header)
        header="#N 4 \n#L Z[cm]  intensityCoh  phaseCoh  intensityEnsemble\n"
        f.write(header)
        
        for i in range(detpoints):
           out = numpy.array((1e2*position3x[i],  fieldIntensity[i], fieldPhase[i], fieldIntensityEA[i]))
           f.write( ("%20.11e "*out.size+"\n") % tuple( out.tolist())  )

    f.close()
    print ("File written to disk: %s"%outfile)
