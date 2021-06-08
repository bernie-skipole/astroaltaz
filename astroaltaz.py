
"""Provides functions to return the Altitude and Azimuth coordinates for an object
   and also the rate of change of Alt and Az of the object.

   The object can be given as a named object such as 'Mars', 'M31', 'Sirius', etc..
   or as equatorial RA, Dec coordinates (which must be provided as decimal degrees float)
   and a timestamp which can be a Python UTC datetime object or Astropy Time object.

   The function "observatory_location" should be edited with the longitude,
   latitude and elevation of the observatory

   NOTE: These functions initiate internet queries, for example to obtain the moon
   or minor planet locations requires a data download

   All functions provided here call upon the astropy/astroquery packages.

   This code is public domain. Feel free to use it as is, or copy and paste what you need.

   Note: astropy uses earth orientation data provided by IERS, it may be useful to have a
   cron job that downloads this weekly to ensure you have a cached up-to date copy. To do this,
   if you have astroplan installed, you can use a short python file containing the lines.

   from astroplan import download_IERS_A

   download_IERS_A()
"""


from datetime import datetime, timedelta

import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, name_resolve, solar_system_ephemeris, get_body, get_sun, get_moon
from astropy.time import Time, TimeDelta
from astroquery.mpc import MPC
from astroquery.exceptions import InvalidQueryError


def observatory_location():
    """Returns the observatory location as an astropy EarthLocation object

    :return: An EarthLocation object
    :rtype: astropy.coordinates.EarthLocation
    """
    # These values should be edited to match the observatory location 
    longitude = -2.1544
    latitude = 53.7111
    elevation = 316
    return EarthLocation.from_geodetic(longitude, latitude, elevation)

# the pressure, temperature, relative_humidity and obswl arguments in the
# following functions are used to calculate the effect of refraction.
# If these arguments are not given, the defaults will ensure no
# refraction calculation is done. All values should be given as floats.

# Pressure is in hPa (equivalent to mbar), if zero, no refraction calculation is made
# Temperature is degrees centigrade
# Relative humidity is a dimensionless quantity between 0 to 1
# obswl is the average wavelength of observations, in micrometers, (1×10−6 metre)

# From the astropy documentation:
# "The refraction model is based on that implemented in ERFA, which is fast but becomes inaccurate for altitudes below about 5 degrees."


def get_named_object(target_name, tstamp, pressure=0.0, temperature=0.0, relative_humidity=0.0, obswl=1):
    """Returns an Astropy SkyCoord object or None if the object cannot be found.
    target_name is a string such as "Mars"
    tstamp is a python datetime or Astropy Time object

    :return: A SkyCoord object
    :rtype: astropy.coordinates.SkyCoord
    """

    if not target_name:
        return

    p = pressure * u.hPa
    t = temperature * u.deg_C
    rh = relative_humidity * u.dimensionless_unscaled
    obswl = obswl * u.micron

    solar_system_ephemeris.set('jpl')
    obsloc = observatory_location()
    if not isinstance(tstamp, Time):
        tstamp = Time(tstamp, format='datetime', scale='utc')
    target_name_lower = target_name.lower()

    aaframe = AltAz(obstime=tstamp, location=obsloc, pressure=p, temperature=t, relative_humidity=rh, obswl=obswl)

    if target_name_lower == "sun":
        target = get_sun(tstamp)
        target_altaz = target.transform_to(aaframe)
        return  target_altaz

    if target_name_lower == "moon":
        target = get_moon(tstamp)
        target_altaz = target.transform_to(aaframe)
        return  target_altaz

    if target_name_lower in ('mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto'):
        target = get_body(target_name_lower, tstamp, obsloc)
        target_altaz = target.transform_to(aaframe)
        return  target_altaz

    # not a planet, see if it is something like M45 or star name
    try:
        target = SkyCoord.from_name(target_name)
    except name_resolve.NameResolveError:
        # failed to find name, maybe a minor planet
        pass
    else:
        target_altaz = target.transform_to(aaframe)
        return  target_altaz

    # minor planet location
    try:
        eph = MPC.get_ephemeris(target_name, location=obsloc, start=tstamp, number=1)
        # eph is a table of a single line, set this into a SkyCoord object
        target = SkyCoord(eph['RA'][0]*u.deg, eph['Dec'][0]*u.deg, obstime = tstamp, location = obsloc, frame='gcrs')
        target_altaz = target.transform_to(aaframe)
    except InvalidQueryError:
        return

    return target_altaz


def get_unnamed_object(target_ra, target_dec, tstamp, pressure=0.0, temperature=0.0, relative_humidity=0.0, obswl=1):
    """Give ra and dec in degrees, return the Astropy SkyCoord object
    target_ra, target_dec are floating point RA, DEC values as decimal degrees 
    tstamp is a datetime or Time object

    :return: A SkyCoord object
    :rtype: astropy.coordinates.SkyCoord
    """

    if (target_ra is None) or (target_dec is None):
        return

    p = pressure * u.hPa
    t = temperature * u.deg_C
    rh = relative_humidity * u.dimensionless_unscaled
    obswl = obswl * u.micron

    solar_system_ephemeris.set('jpl')
    obsloc = observatory_location()
    if not isinstance(tstamp, Time):
        tstamp = Time(tstamp, format='datetime', scale='utc')

    target = SkyCoord(target_ra*u.deg, target_dec*u.deg, frame='icrs')
    return target.transform_to(AltAz(obstime=tstamp, location=obsloc, pressure=p, temperature=t, relative_humidity=rh, obswl=obswl))


def get_named_alt_az(target_name, tstamp, pressure=0.0, temperature=0.0, relative_humidity=0.0, obswl=1):
    """Return alt, az, as decimal degrees, given a target name
    and time stamp (python utc datetime or astropy Time object) 
    If nothing found, returns None

    :return: Altitude, Azimuth, as decimal degrees
    :rtype: Tuple of two Floats
    """

    target_altaz = get_named_object(target_name, tstamp, pressure, temperature, relative_humidity, obswl)
    if target_altaz is None:
        return
    return target_altaz.alt.degree, target_altaz.az.degree


def get_unnamed_alt_az(target_ra, target_dec, tstamp, pressure=0.0, temperature=0.0, relative_humidity=0.0, obswl=1):
    """Return alt, az, as decimal degrees, given J2000 coordinates target_ra, target_dec which should be
    decimal degrees and time stamp which can be either a python utc datetime or astropy Time object.

    :return: Altitude, Azimuth, as decimal degrees
    :rtype: Tuple of two Floats
    """
    target_altaz = get_unnamed_object(target_ra, target_dec, tstamp, pressure, temperature, relative_humidity, obswl)
    if target_altaz is None:
        return
    return target_altaz.alt.degree, target_altaz.az.degree


def get_named_alt_az_rates(target_name, tstamp, pressure=0.0, temperature=0.0, relative_humidity=0.0, obswl=1):
    """Return alt, az, alt_rate, az_rate, given a target name
    and time stamp which can be a python utc datetime or astropy Time object
    alt,az returned will be decimal degrees, the rates returned are degrees per second
    If nothing found, returns None

    :return: Alt, Az, Alt deg per sec, Az deg per sec
    :rtype: Tuple of four Floats
    """

    if not isinstance(tstamp, Time):
        tstamp = Time(tstamp, format='datetime', scale='utc')
    td = TimeDelta(10, format='sec')
    coords = get_named_alt_az(target_name, tstamp, pressure, temperature, relative_humidity, obswl)
    if coords is None:
        return
    # get coordinates 10 seconds in the past
    # and 10 seconds in the future
    # and set the rates as degrees changed / 20
    coords_minus = get_named_alt_az(target_name, tstamp-td, pressure, temperature, relative_humidity, obswl)
    coords_plus = get_named_alt_az(target_name, tstamp+td, pressure, temperature, relative_humidity, obswl)
    # altitude rate
    alt_rate = (coords_plus[0] - coords_minus[0])/20.0
    # azimuth rate
    # handle case where azimuth crosses the 360 - 0 border
    if coords_minus[1]>270 and coords_plus[1]<90:                   # example  cp = 2, cm = 358
        az_rate = (360 + coords_plus[1] - coords_minus[1])/20.0     # rate = (360 + 2 - 358)/20 = 4/20
    elif coords_minus[1]<90 and coords_plus[1]>270:                 # example cp = 358, cm = 2
        az_rate = (coords_plus[1] - 360 - coords_minus[1])/20.0     # rate = (358 - 360 - 2)/20 = -4/20
    else:
        az_rate = (coords_plus[1] - coords_minus[1])/20.0
    return coords[0], coords[1], alt_rate, az_rate


def get_unnamed_alt_az_rates(target_ra, target_dec, tstamp, pressure=0.0, temperature=0.0, relative_humidity=0.0, obswl=1):
    """Return alt, az, alt_rate, az_rate, given J2000 coordinates target_ra and target_dec which should be provided as
    floats of decimal degrees, and time stamp (python utc datetime or astropy Time object)
    alt,az returned will be decimal degrees, the rates returned are degrees per second

    :return: Alt, Az, Alt deg per sec, Az deg per sec
    :rtype: Tuple of four Floats
    """
    if not isinstance(tstamp, Time):
        tstamp = Time(tstamp, format='datetime', scale='utc')
    td = TimeDelta(10, format='sec')
    coords = get_unnamed_alt_az(target_ra, target_dec, tstamp, pressure, temperature, relative_humidity, obswl)
    if coords is None:
        return
    # get coordinates 10 seconds in the past
    # and 10 seconds in the future
    # and set the rates as degrees changed / 20
    coords_minus = get_unnamed_alt_az(target_ra, target_dec, tstamp-td, pressure, temperature, relative_humidity, obswl)
    coords_plus = get_unnamed_alt_az(target_ra, target_dec, tstamp+td, pressure, temperature, relative_humidity, obswl)
    # altitude rate
    alt_rate = (coords_plus[0] - coords_minus[0])/20.0
    # azimuth rate
    # handle case where azimuth crosses the 360 - 0 border
    if coords_minus[1]>270 and coords_plus[1]<90:                   # example  cp = 2, cm = 358
        az_rate = (360 + coords_plus[1] - coords_minus[1])/20.0     # rate = (360 + 2 - 358)/20 = 4/20
    elif coords_minus[1]<90 and coords_plus[1]>270:                 # example cp = 358, cm = 2
        az_rate = (coords_plus[1] - 360 - coords_minus[1])/20.0     # rate = (358 - 360 - 2)/20 = -4/20
    else:
        az_rate = (coords_plus[1] - coords_minus[1])/20.0
    return coords[0], coords[1], alt_rate, az_rate


if __name__ == "__main__":
    # Examples - for various objects, print ten lines of ra,dec,ra rate, dec rate for each at 15 minute intervals
    tstamp = datetime.utcnow()
    td = timedelta(minutes=15)
    # test with Mizar
    print("Mizar")
    for n in range(10):
        coords = get_named_alt_az_rates("Mizar", tstamp)
        print(coords)
        tstamp += td
    # test with arbitrary RA DEC point
    tstamp = datetime.utcnow()
    print("RA 80.0 deg, DEC 60.0 deg")
    for n in range(10):
        coords = get_unnamed_alt_az_rates(80.0, 60.0, tstamp)
        print(coords)
        tstamp += td
    # test with the Sun
    tstamp = datetime.utcnow()
    print("Sun")
    for n in range(10):
        coords = get_named_alt_az_rates("Sun", tstamp)
        print(coords)
        tstamp += td
    # test with planet
    tstamp = datetime.utcnow()
    print("Mars")
    for n in range(10):
        coords = get_named_alt_az_rates("Mars", tstamp)
        print(coords)
        tstamp += td
    # test with minor planet Ceres
    tstamp = datetime.utcnow()
    print("Ceres")
    for n in range(10):
        coords = get_named_alt_az_rates("Ceres", tstamp)
        print(coords)
        tstamp += td
    # test with the Moon
    tstamp = datetime.utcnow()
    print("Moon")
    for n in range(10):
        coords = get_named_alt_az_rates("Moon", tstamp)
        print(coords)
        tstamp += td
    


