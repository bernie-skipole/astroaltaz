The astroaltaz module
=====================

.. automodule:: astroaltaz

The functions
^^^^^^^^^^^^^

.. autofunction:: astroaltaz.observatory_location

This function "observatory_location" should be edited with the longitude, latitude and elevation of the observatory.


In the following functions the pressure, temperature, relative_humidity and obswl arguments are used to calculate the effect of refraction.

If these arguments are not given, the defaults will ensure no refraction calculation is done

Pressure is in hPa (equivalent to mbar), if zero, no refraction calculation is made

Temperature is degrees centigrade

Relative humidity is a dimensionless quantity between 0 to 1

obswl is the average wavelength of observations, in micrometers, (1×10−6 metre)

From the astropy documentation:

"The refraction model is based on that implemented in ERFA, which is fast but becomes inaccurate for altitudes below about 5 degrees."


.. autofunction:: astroaltaz.get_named_object

.. autofunction:: astroaltaz.get_unnamed_object

.. autofunction:: astroaltaz.get_named_alt_az

.. autofunction:: astroaltaz.get_unnamed_alt_az

.. autofunction:: astroaltaz.get_named_alt_az_rates

.. autofunction:: astroaltaz.get_unnamed_alt_az_rates



