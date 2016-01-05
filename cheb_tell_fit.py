#! /usr/bin/env python

"""
Sam Vaughan: 05/012016

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
"""
from __future__ import print_function

import numpy as np
import numpy.ma as ma
from numpy.polynomial.chebyshev import chebfit,chebval
import matplotlib.pyplot as plt
from astropy.io import fits
import argparse

from ast import literal_eval


parser = argparse.ArgumentParser(description="Fit Chebyshev polynomials to parts of a telluric spectrum away from sky features. Called as final_tell.py <inputfilename> <outputfilename>")
parser.add_argument("inputfilename")
parser.add_argument("outputfilename")
parser.add_argument("--show", help="Display a plot of the original spectrum, a poly-fit smoothed version and a non smoothed version rather than saving", action='store_true')
parser.add_argument("--order", help="Change the order of the Chebyshev polynomial fitting. Default is 10", type=int)
parser.add_argument("--pixels", help="Change the number of pixels to blend at the start and end of each polynomial/spectra join. Default is 10", type=int)
parser.add_argument("--tellrange", nargs=1, type=str, help="A file with the wavelength ranges of telluric features you want to fit, in the form [lamda1, lamda2]\n[lamda3, lamda4]\n etc]")
args = parser.parse_args()
filename=args.inputfilename
outfile=args.outputfilename


#Defaults
show=False
order=10
pixels=10

#Examples of telluric features range
wavelength_ranges=[[0.638, 0.66], [0.686, 0.694], [0.705, 0.74], [0.759, 0.775], [0.810, 0.843], [0.891, 0.984]]
if args.show:
    show=True

if args.order:
    order=args.order

if args.pixels:
    pixels=args.pixels

if args.tellrange:

    wavelength_ranges=[]
    tell_filename=args.tellrange[0]
    print("Using telluric features from file {}".format(tell_filename))
    with open(tell_filename, "r") as f:
        for line in f:
            wavelength_ranges.append(literal_eval(line.strip()))


print("Telluric Features are: {}".format(wavelength_ranges))

hdulist=fits.open(filename)
spec=hdulist[0].data

crpix=hdulist[0].header["CRPIX1"]
crval=hdulist[0].header["CRVAL1"]
cdelt=hdulist[0].header["CDELT1"]

lamdas=np.array([(l-crpix)*cdelt +crval for l in range(len(spec))])




print("Polynomial Fitting order is {}, Blending is over {} pixels".format(order, pixels))


def start_blend(start, pixels, fitted_spec, spectrum):
    """Blend the polynomial and spectrum at a spectrum-polynomial join"""
    for i in range(pixels):
        fitted_spec[start+i]=(1-0.1*i)*spectrum[start+i]+0.1*i*fitted_spec[start+i]
    return fitted_spec[start:start+pixels]

def end_blend(start, pixels, fitted_spec, spectrum):
    """Blend the polynomial and spectrum at a polynomial-spectrum join"""
    for i in range(pixels):
        fitted_spec[start-i]=0.1*i*fitted_spec[start-i]+(1-0.1*i)*spectrum[start-i]
    return fitted_spec[start-pixels:start]




#List of indices correspoding to telluric features
index_ranges=[[np.where(lamdas==val[0]), np.where(lamdas==val[1])] for val in wavelength_ranges]


#Copies of the input spectra.
#We'll use one to plot the polynomials with no smoothing, one with
fit_spec=np.copy(spec)
fit_spec_no_smooth=np.copy(spec)


start=0
for i, index in enumerate(index_ranges):

        #Where the telluric features are
        t_start=int(index[0][0])
        t_end=int(index[1][0])



        #Get the indices away from the telluric features
        end=t_start

        #Rescale lamda array to be between -1 and 1
        l=np.array([2*(w-lamdas[start:end].min())/(lamdas[start:end].max()-lamdas[start:end].min())-1 for w in lamdas[start:end]])

        polyfit=chebfit(l, spec[start:end], order)
        polynomial=chebval(l, polyfit)

        print("Fitting Polynomials from lamda={}A to lamda={}A".format(lamdas[start], lamdas[end]))
        fit_spec[start:end]=polynomial
        fit_spec_no_smooth[start:end]=polynomial
        fit_spec[start:start+pixels]=start_blend(start, pixels, fit_spec, spec)
        fit_spec[end-pixels:end]=end_blend(end, pixels, fit_spec, spec)


        #Update the start value for the next section of the spectrum to fit
        start=t_end

#Do the last section of the fit
s=index_ranges[-1][1][0]
#Rescale the wavelength array to be between -1 and 1.
l=[2*(w-lamdas[s:].min())/(lamdas.max()-lamdas.min())-1 for w in lamdas[s:]]
polyfit=chebfit(l, spec[s:], order)
polynomial=chebval(l, polyfit)
fit_spec[s:]=polynomial
fit_spec_no_smooth[s:]=polynomial

fit_spec[start:start+pixels]=start_blend(start, pixels, fit_spec, spec)

if show==False:
    print("Writing to {}".format(outfile))
    hdulist[0].data=fit_spec

    hdulist.writeto(outfile)


if show==True:
    plt.plot(lamdas, fit_spec_no_smooth, c="r")
    plt.plot(lamdas, fit_spec, c="b")

    plt.plot(lamdas, spec, c="g")
    #.plot(lamdas, spec)
    plt.show()
