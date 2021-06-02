# astroaltaz
Python functions to return the Altitude and Azimuth coordinates for an object and also the rate of change of Alt and Az of the object.

The object can be given as a named object such as 'Mars', 'M31', 'Sirius', etc.. or as equatorial RA, Dec coordinates (which must be provided as decimal degrees float)  and a timestamp which can be a Python UTC datetime object or Astropy Time object.

The function "observatory_location" should be edited with the longitude, latitude and elevation of the observatory.

NOTE: These functions initiate internet queries, for example to obtain the moon or minor planet locations requires a data download.

All functions provided here use Python3 and call upon the astropy and astroquery packages.

The code in this file is public domain. Feel free to use it as is, or copy and paste what you need. The file is very short, since all the actual work is done by astropy/astroquery, please look at the source which is well commented.


Note: astropy uses earth orientation data provided by IERS, it may be useful to have a cron job that downloads this weekly to ensure you have a cached up-to date copy. To do this, if you have astroplan installed, you can use a short python file containing the lines.


from astroplan import download_IERS_A

download_IERS_A()



