# astroaltaz
Python functions to return the Altitude and Azimuth coordinates for an object and also the rate of change of Alt and Az of the object.

The object can be given as a named object such as 'Mars', 'M31', 'Sirius', etc.. or as equatorial RA, Dec coordinates (which must be provided as decimal degrees float)  and a timestamp which can be a Python UTC datetime object or Astropy Time object.

The function "observatory_location" should be edited with the longitude, latitude and elevation of the observatory.

NOTE: These functions initiate internet queries, for example to obtain the moon or minor planet locations requires a data download.

All functions provided here call upon the astropy/astroquery packages.

This code is public domain. Feel free to use it as is, or copy and paste what you need.
