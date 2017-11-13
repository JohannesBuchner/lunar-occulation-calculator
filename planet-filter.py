"""
Detailed filter for lunar occultations from a candidate list.

For each candidate, it does the following:

- search for the closest angular separation (between target and moon center).
- from there, search for the entrance and exit of the occultation.
- measure how long a target circle of 1'' remains partially occulted

The latter quantity tells you whether the occultation is more head-on or grazing.
A grazing lunar occultation is better because the occultation is slower.

You need to run lunar-filter-rough.py first.

Johannes Buchner (C) 2017
"""

import sys
import numpy
import matplotlib.pyplot as plt
import astropy
import astropy.coordinates
import astropy.units as u
import scipy.optimize
import joblib

mem = joblib.Memory('.', verbose=False)

lon_deg = 107.
lat_deg = 34.
height = 2000. * u.m
lon = astropy.coordinates.Longitude(lon_deg * u.deg)
lat = astropy.coordinates.Latitude(lat_deg * u.deg)
loc = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)

@mem.cache
def compute_body_separation(time, target, body):
	planet = astropy.coordinates.get_body(body=body, time=time, location=loc)
	return planet.separation(target)

# get a list of nearby AGN
started = False
step = astropy.time.TimeDelta(u.hour * 100000, format='sec')

for targeti, line in enumerate(open('candidates_planets')):
	ra, dec, approxtime, body = line.strip().split('\t')
	ra, dec, approxtime = float(ra), float(dec), float(approxtime)
	target = astropy.coordinates.SkyCoord(ra, dec, unit=u.deg, frame=astropy.coordinates.FK5, equinox='J2000')
	
	if body in ('saturn', 'mars'):
		test_radius = 30 * u.arcsec
	elif body in ('uranus', 'neptune'):
		test_radius = 5 * u.arcsec
	elif body in ('venus', 'jupiter'):
		test_radius = 60 * u.arcsec
	else:
		assert False, (body, 'angular size unknown')
	starttime = astropy.time.Time(approxtime, format='jd').utc
	
	def minfunc((i,)):
		time = starttime + i*step
		separation = compute_body_separation(time, target, body)
		return separation.arcsec
	
	iopt, = scipy.optimize.fmin_cobyla(minfunc, [0.], [], disp=0)
	time = starttime + iopt*step
	separation = compute_body_separation(time, target, body)
	print 'closest to %s: %.3f arcsec [need %s]' % (body, separation.arcsec, test_radius), time
	if separation > test_radius:
		print '     --> skipping.'
		continue
	

