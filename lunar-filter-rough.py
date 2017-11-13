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
# TODO: vary
lon = astropy.coordinates.Longitude(107 * u.deg)
lat = astropy.coordinates.Latitude(34 * u.deg)
height = 2124. * u.m
loc = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)

starttime = astropy.time.Time.now()

# make a time step such that the moon advances a fraction of its angular diameter
fraction = 1./2.
moon_radius = 31. * u.arcmin / 2.
test_radius = 10 * moon_radius
moon_speed = (27. * 24 * u.hour) / (360 * u.deg) # very approximately
moon_step = (moon_speed * test_radius.to(u.deg)).to(u.second)
step = astropy.time.TimeDelta(moon_step * fraction, format='sec')
print 'time step:', step
candidates = set()
candidates_file = open('candidates', 'w')

for i in range(100000):
	time = starttime + i*step
	moon = astropy.coordinates.get_moon(time, location=loc)
	separation = moon.separation(targets)
	sys.stderr.write('%d ... %s %.3f arcmin  \r' % (i, time, separation.min().arcmin))
	# now we are covering something
	occulted = separation < test_radius # roughly, to get candidate list
	if not occulted.any():
		continue
	
	for i in numpy.where(occulted)[0]:
		if i not in candidates:
			candidates.add(i)
			candidates_file.write('%.10f\t%.10f\t%.1f\n' % (ras[i], decs[i], time.jd))
			candidates_file.flush()
			print 'new candidate found:', i, names[i], moon.separation(targets)[i], time



