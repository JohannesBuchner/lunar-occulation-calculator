"""
Detailed filter for lunar occultations from a candidate list.

For each candidate, it does the following:

- search for the closest angular separation (between target and moon center).
- from there, search for the entrance and exit of the occultation.
- measure how long a target circle of 1'' remains partially occulted

The latter quantity tells you whether the occultation is more head-on or grazing.
A grazing lunar occultation is better because the occultation is slower.

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

@mem.cache
def compute_moon_separation(lon_deg, lat_deg, height, time, target):
	lon = astropy.coordinates.Longitude(lon_deg * u.deg)
	lat = astropy.coordinates.Latitude(lat_deg * u.deg)
	loc = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)

	moon = astropy.coordinates.get_moon(time, location=loc)
	return moon.separation(target)

# get a list of nearby AGN
started = False
ras, decs, approxtimes = numpy.loadtxt('candidates', ndmin=2).transpose()
# make a time step such that the moon advances a fraction of its angular diameter
fraction = 1./10.
moon_radius = 31. * u.arcmin / 2.
test_radius = 1 * moon_radius
moon_speed = (27. * 24 * u.hour) / (360 * u.deg) # very approximately
moon_step = (moon_speed * test_radius.to(u.deg)).to(u.second)
step = astropy.time.TimeDelta(moon_step * fraction, format='sec')
print 'time step:', step
print 'moon speed:', 1/(moon_speed.to(u.second / u.arcsec))

for targeti, (ra, dec, approxtime) in enumerate(zip(ras, decs, approxtimes)):
	target = astropy.coordinates.SkyCoord(ra, dec, unit=u.deg, frame=astropy.coordinates.FK5, equinox='J2000')

	starttime = astropy.time.Time(approxtime, format='jd')
	
	height = 2000. * u.m
	def minfunc((i, lon_deg, lat_deg)):
		# observer:
		if not 0 < lon_deg < 360:
			return 1e10
		if not -90 < lat_deg < 90:
			return 1e10
		#lon = astropy.coordinates.Longitude(lon_deg * u.deg)
		#lat = astropy.coordinates.Latitude(lat_deg * u.deg)
		#loc = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)

		time = starttime + i*step*10000
		separation = compute_moon_separation(lon_deg, lat_deg, height, time, target)
		#moon = astropy.coordinates.get_moon(time, location=loc)
		#separation = moon.separation(target)
		#print 'separation: %.2f' % separation.arcmin, i, lon_deg, lat_deg
		return separation.arcmin
	iopt, lon_deg, lat_deg = scipy.optimize.fmin_cobyla(minfunc, [0., 107., 34.], [], disp=0)
	lon = astropy.coordinates.Longitude(lon_deg * u.deg)
	lat = astropy.coordinates.Latitude(lat_deg * u.deg)
	time = starttime + iopt*step*10000
	loc = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)
	moon = astropy.coordinates.get_moon(time, location=loc)
	separation = moon.separation(target)
	print 'best separation: %.3f arcmin' % separation.arcmin, time
	if separation > moon_radius:
		print ' --> useless, not occulted. skipping.'
		continue
	
	#continue
	
	# now the question is: if we have a disk of size moon_radius
	# and a target disk of size 1 arcsec
	# how long is the target disk partially occulted
	# here we should also optimize the location
	
	# we are looking for a time where the distance between moon center and target is almost the moon radius
	# this should happen both before and after the current best time, which is where moon center ~= target
	# we want to compute the moon angular velocity, i.e. by what fraction the occultation
	# of the target disk changed. This should be minimized.

	for sign in -1, +1:
	
		def minfunc((i, lon_deg, lat_deg)):
			# observer:
			if i * sign > 0: 
				return 1e20 + sign * 1e20
			if not 0 < lon_deg < 360:
				return 1e20
			if not -90 < lat_deg < 90:
				return 1e20
			#lon = astropy.coordinates.Longitude(lon_deg * u.deg)
			#lat = astropy.coordinates.Latitude(lat_deg * u.deg)
			#loc = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)

			time = starttime + (i+iopt)*step*10000
			#moon = astropy.coordinates.get_moon(time, location=loc)
			#separation = moon.separation(target)
			separation = compute_moon_separation(lon_deg, lat_deg, height, time, target)
			# now compute the separation between target disk and moon
			dist_beyond_moon = separation - moon_radius
			# coverage is zero if |dist_beyond_moon| < 1arcsec
			quality = 1e10 + 1e10 * (dist_beyond_moon.arcsec)**2
			if numpy.abs(dist_beyond_moon.arcsec) > 1.:
				print '%.2f    free range: %.2f arcsec (should be <1arcsec)' % (quality, dist_beyond_moon.arcsec), i, lon_deg, lat_deg
				return quality
		
			# now we are partially occulted by the moon.
			# technique:
			# advance 1000s and measure new separation. We want this to be slow.
			time_later = time + 1000 * u.second
			#moon_later = astropy.coordinates.get_moon(time_later, location=loc)
			#separation_later = moon_later.separation(target)
			separation_later = compute_moon_separation(lon_deg, lat_deg, height, time_later, target)
			dist_beyond_moon_later = separation_later - moon_radius
			# we want to minimize the difference in moon separation
			quality = (dist_beyond_moon.arcmin)**2 + (dist_beyond_moon - dist_beyond_moon_later).arcmin**2
			print '%.2f  free range: %.2f .. %.2f arcsec' % (quality, dist_beyond_moon.arcsec, dist_beyond_moon_later.arcsec), i, lon_deg, lat_deg
			return quality
	
		iopt2, lon_deg2, lat_deg2 = scipy.optimize.fmin(minfunc, [0, lon_deg, lat_deg])
		
		lon = astropy.coordinates.Longitude(lon_deg2 * u.deg)
		lat = astropy.coordinates.Latitude(lat_deg2 * u.deg)
		loc = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)

		time = starttime + (iopt2+iopt)*step*10000
		
		print 'computing moon grid ...'
		#times = time + numpy.linspace(-1000, 1000, 1000) * u.second
		#separation = [compute_moon_separation(lon_deg, lat_deg, height, time_later, target) for time_later in times]
		times = []
		separations = []
		for i in range(0, 1000, 10):
			time_later = time + i * u.second
			separation1 = compute_moon_separation(lon_deg, lat_deg, height, time_later, target)
			times.append(i)
			separations.append(separation1.arcsec)
			time_later = time - i * u.second
			separation2 = compute_moon_separation(lon_deg, lat_deg, height, time_later, target)
			times.insert(0, -i)
			separations.insert(0, separation2.arcsec)
			print i, (separation1 - moon_radius).arcsec, (separation2 - moon_radius).arcsec
			if numpy.abs((separation1 - moon_radius).arcsec) > 10 and numpy.abs((separation2 - moon_radius).arcsec) > 10:
				break
		
		x = time + numpy.array(times) * u.second
		dt = numpy.array(times)
		y = numpy.array(separations) * u.arcsec - moon_radius
		print 'storing solutions ...'
		filename = 'occultation_%d%s' % (targeti, 'p' if sign > 0 else 'm')
		numpy.savetxt(filename + '.txt', numpy.transpose([x.jd, separations]))
		print 'plotting ...', moon_radius
		plt.plot(dt, y.value, color='k')
		plt.plot(dt, (dt / moon_speed.to(u.second / u.arcsec)).value, ':', color='gray')
		plt.ylim(y.value.min(), y.value.max())
		plt.savefig(filename + '.pdf', bbox_inches='tight')
		plt.close()
	

