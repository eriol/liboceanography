Function reference
==================

``liboceanography`` define the following functions:

#. adiabatic_temperature_gradient
#. atg
#. conductivity
#. cpsw
#. depth
#. freezing_point
#. potential_temperature
#. salinity
#. sound_speed
#. specific_heat
#. specific_volume_anomaly
#. svan
#. theta


.. _ref_adiabatic_temperature_gradient:

adiabatic_temperature_gradient
------------------------------

Compute the adiabatic temperature gradient.

.. code-block:: c

    double adiabatic_temperature_gradient(double salinity, double temperature,
                                          double pressure)

Units::

    salinity -- PSS-78
    temperature -- degrees Celsius
    pressure  -- decibars

Returns adiabatic temperature gradient in °C/decibars.

atg
---

Alias for :ref:`ref_adiabatic_temperature_gradient`.

conductivity
------------

Convert salinity to conductivity ratio.

.. code-block:: c

    double conductivity(double salinity, double temperature, double pressure)

Units::

    salinity -- PSS-78
    temperature -- degrees Celsius
    pressure  -- decibars

cpsw
----

Alias for :ref:`ref_specific_heat`.

depth
-----

Compute depth from pressure using Saunders and Fofonoff's method.

.. code-block:: c

    double depth(double pressure, double latitude)

Units::

    pressure  -- decibars
    latitude -- degrees

Returns depth in meters.

freezing_point
--------------

Compute the freezing point of seawater.

.. code-block:: c

    double freezing_point(double salinity, double pressure)

Units::

    salinity -- PSS-78
    pressure  -- decibars

Returns freezing point in degrees Celsius.

.. _ref_potential_temperature:

potential_temperature
---------------------

Compute the local potential temperature at reference pressure.

.. code-block:: c

    double potential_temperature(double salinity, double temperature,
                                 double pressure, double reference_pressure)

Units::

    salinity -- PSS-78
    temperature -- degrees Celsius
    pressure  -- decibars
    reference_pressure  -- decibars

Returns local potential temperature in degrees Celsius.

salinity
--------

Convert conductivity ratio to salinity.

.. code-block:: c

    double salinity(double conductivity, double temperature, double pressure)

Units::

    temperature -- degrees Celsius
    pressure  -- decibars

Returns salinity in PSS-78.

sound_speed
-----------

Compute the speed of sound in seawater by Chen and Millero.

.. code-block:: c

    double sound_speed(double salinity, double temperature, double pressure)

Units::

    salinity -- PSS-78
    temperature -- degrees Celsius
    pressure  -- decibars

Returns sound speed in meters/second.

.. _ref_specific_heat:

specific_heat
-------------

Compute the specific heat of seawater.

.. code-block:: c

    double specific_heat(double salinity, double temperature, double pressure)

Units::

    salinity -- PSS-78
    temperature -- degrees Celsius
    pressure  -- decibars

Returns specific heat in J/(Kg °C).

.. _ref_specific_volume_anomaly:

specific_volume_anomaly
-----------------------

Compute specific volume anomaly (steric anomaly).

.. code-block:: c

    double specific_volume_anomaly(double salinity, double temperature,
                                double pressure, double *sigma)

Units::

    salinity -- PSS-78
    temperature -- degrees Celsius
    pressure  -- decibars
    sigma (density anomaly) -- Kg/m^3

Returns specific volume anomaly as 1.0e-8 m^3/Kg.

svan
----

Alias for :ref:`ref_specific_volume_anomaly`.

theta
-----

Alias for :ref:`ref_potential_temperature`.
