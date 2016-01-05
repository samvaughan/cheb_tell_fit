## cheb_tell_fit.py
#Fit Chebsyshev polynomials to sections of an A0V telluric spectrum away from sky features, for a less noisy telluric correction.


A script to fit chebyshev polynomials to the parts of a telluric spectrum (already divided by an A0 star spec using matchtelluric.pro)
away from the telluric features, and keep the observed spectrum everywhere else. The polynomials are blended in with the raw spectrum for
(by default) 10 pixels at either join.

Called as cheb_tell_fit.py <inputfilename> <outputfilename>

Optional arguments:

--show: Rather than saving the output, display a plot of the original spectrum, a poly-fit smoothed version and a non smoothed version
--order: Change the order of the Chebyshev polynomial fitting. Default is 10
--pixels: Change the number of pixels to blend at the start and end of each polynomial/spectra join. Default is 10
--tellrange: A text file with a list of wavelength ranges, each containing a section of the spectrum you want to keep. Polynomials will be fit in between these ranges. They should be in the
format [lamda1, lamda2]\n[lamda3,lamda4]\n...etc. The square brackets are important as they're read into python using literal_eval (i.e read to be a two element list straight away).


Type python cheb_tell_fit -h for info on the options

Sam Vaughan
