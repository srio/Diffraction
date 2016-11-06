"""


correlation_lengths.py: calculates


"""

import numpy


__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2016"

def std0(x):
    mu = 0.0 # x.mean()
    return numpy.sqrt(numpy.mean(abs(x - mu)**2))

def coherence_length_longitudinal(lambda1,lambda2,lambda0=None):
    delta_lambda = numpy.abs(lambda1 - lambda2)
    if lambda0 == None:
        lambda0 = 0.5 * numpy.abs(lambda1 + lambda2)
    return lambda0**2 / (2 * delta_lambda)

def coherence_length_transverse(lambda0,distance,separation_at_source):
    return lambda0 * distance / (2 * numpy.abs(separation_at_source))

def coherence_length_histogram(lambda0,distance,sigma,mu=0.0,npoints=1000,do_plot=True):
    #source
    s1 = numpy.random.normal(mu, sigma, npoints)
    s2 = numpy.random.normal(mu, sigma, npoints)

    # s1 = (numpy.random.random(npoints) - 0.5 ) * sigma * 2.35
    # s2 = (numpy.random.random(npoints) - 0.5 ) * sigma * 2.35
    # separation_at_source
    x1 = numpy.outer(s1,numpy.ones(s2.size))
    x2 = numpy.outer(numpy.ones(s1.size),s2)
    dd = numpy.abs(x1-x2).flatten()
    separation_at_source = numpy.unique(dd)[1:-1]

    cl = coherence_length_transverse(lambda0,distance,separation_at_source)

    if do_plot:
        import matplotlib.pyplot as plt

        count, bins, ignored = plt.hist(s1, 100, normed=True)
        plt.plot(bins, 1/(sigma * numpy.sqrt(2 * numpy.pi)) * numpy.exp( - (bins - mu)**2 / (2 * sigma**2) ),
                 linewidth=2, color='r')
        plt.show()


        # https://arxiv.org/pdf/1508.02238v1.pdf
        count, bins, ignored = plt.hist(1e6*separation_at_source, 100, normed=True)
        sigma_estimated = sigma * numpy.sqrt(2)
        print("StDev of source separation = %f um; estimating = %f um"%
              (1e6*std0(separation_at_source),1e6*sigma_estimated))
        plt.plot(bins, 2/(1e6*sigma_estimated * numpy.sqrt(2 * numpy.pi)) * numpy.exp( - (bins - 0.0)**2 / (2 * (1e6*sigma_estimated)**2) ),
                 linewidth=2, color='r')
        plt.show()


        #
        count, bins, ignored = plt.hist(1e6*cl, 500, normed=True, range = [0,10*1e6*coherence_length_transverse(lambda0,distance,2.35*sigma)])
        plt.show()

    return bins[count.argmax()]*1e-6




if __name__ == "__main__":
    lambda0 = 1e-10
    distance = 30.0

    delta_lambda = lambda0 * 1e-4
    lambda1 = lambda0 - 0.5 * delta_lambda
    lambda2 = lambda0 + 0.5 * delta_lambda

    sigma_v = 3.5e-6
    sigma_lb = 37.4e-6
    sigma_hb = 387.8e-6
    sigma_ebs = 1e-6 # 27.2e-6
    print("Wavelength = %f A"%(lambda0*1e10))

    print("Correlation length (longitudinal): %f um"%(1e6*coherence_length_longitudinal(lambda1,lambda2)))
    print("Correlation length V     (transversal): %f um"%(1e6*coherence_length_transverse(lambda0,distance,2.35*sigma_v)))
    print("Correlation length H Hb  (transversal): %f um"%(1e6*coherence_length_transverse(lambda0,distance,2.35*sigma_hb)))
    print("Correlation length H Lb  (transversal): %f um"%(1e6*coherence_length_transverse(lambda0,distance,2.35*sigma_lb)))
    print("Correlation length H EBS (transversal): %f um"%(1e6*coherence_length_transverse(lambda0,distance,2.35*sigma_ebs)))


    # print("HISTOGRAM Correlation length (longitudinal): %f um"%(1e6*coherence_length_longitudinal(lambda1,lambda2)))
    # print("HISTOGRAM Correlation length V     (transversal): %f um"%(1e6*coherence_length_histogram(lambda0,distance,sigma_v)))
    # print("HISTOGRAM Correlation length H Hb  (transversal): %f um"%(1e6*coherence_length_histogram(lambda0,distance,sigma_hb)))
    # print("HISTOGRAM Correlation length H Lb  (transversal): %f um"%(1e6*coherence_length_histogram(lambda0,distance,sigma_lb)))
    print("HISTOGRAM Correlation length H EBS (transversal): %f um"%(1e6*coherence_length_histogram(lambda0,distance,sigma_ebs)))

