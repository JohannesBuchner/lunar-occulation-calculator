"""
Rough filter for lunar occultations with a target list.

Johannes Buchner (C) 2017
"""


import sys
import numpy
import astropy
import astropy.coordinates
import astropy.units as u

# get a list of nearby AGN
started = False
ras = []
decs = []
names = []
for line in open('swiftbatcoords.tsv'):
	if line.startswith('----'):
		started = True
		continue
	if not started: continue
	parts = line.split('\t')
	ras.append(float(parts[0]))
	decs.append(float(parts[1]))
	names.append(parts[3])

targets = astropy.coordinates.SkyCoord(ras, decs, frame=astropy.coordinates.FK5, unit=u.deg, equinox='J2000')

# observer:
lon = astropy.coordinates.Longitude(107 * u.deg)
lat = astropy.coordinates.Latitude(34 * u.deg)
height = 2124. * u.m
loc = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)

starttime = astropy.time.Time.now()
test_radius = 1. * u.deg
timestep = 30. * u.hour

candidates = set()
candidates_file = open('candidates_planets', 'w')

for i in range(100000):
	time = starttime + i*timestep
	bodies = ['venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
	for body in bodies:
		planet = astropy.coordinates.get_body(body=body, time=time, location=loc)
		separation = planet.separation(targets)
		sys.stderr.write('%d ... %s %.3f arcmin  \r' % (i, time, separation.min().arcmin))
		# now we are covering something
		occulted = separation < test_radius # roughly, to get candidate list
		if not occulted.any():
			continue
	
		for i in numpy.where(occulted)[0]:
			if i not in candidates:
				candidates.add(i)
				candidates_file.write('%.10f\t%.10f\t%.1f\t%s\n' % (ras[i], decs[i], time.jd, body))
				candidates_file.flush()
				print 'new candidate found:', i, names[i], separation[i], time, body



