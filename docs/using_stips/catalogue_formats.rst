STIPS Catalog Formats
=====================
.. note::

    If you use the STIPS Scene Module to create your observation scenes, then
    created catalogs will already be in supported formats (specifically the
    Phoenix format for stars and the BC95 format for galaxies). Alternately, if
    you create new STIPS catalogs by editing existing catalogs, then, as long
    as you keep to the existing format, any created catalogs should work.

General Formatting
------------------

In general, STIPS catalogs may be in one of two formats, IPAC_ text or FITS
BinTable_. Both are accessed through the `Astropy Table
<https://docs.astropy.org/en/stable/table/>`_ API. In the case of FITS tables,
STIPS will internally divide tables into chunks of 100,000 entries, putting each
chunk into its own extension.

In IPAC_ tables, metadata is presented at the start of the table, in lines with
the format::

    \key=value

In FITS tables, metadata is included as header keywords.

STIPS Table Formats
-------------------

* :ref:`phoenix-tables`

* :ref:`realtime-phoenix-tables`

* :ref:`pandeia-tables`

* :ref:`bc95-catalog`

* :ref:`internal-catalog`

*  :ref:`mixed-catalog`

* :ref:`multifilter-catalog`

* :ref:`generic-catalog`

* Output

.. _phoenix-tables:

Phoenix Tables
--------------

Phoenix tables are tables of point sources that can be generated as models using
the `Phoenix <http://phoenix.ens-lyon.fr>`_ star simulator. STIPS uses
`stsynphot <https://stsynphot.readthedocs.io/en/latest/>`_ to generate its
Phoenix spectra, either using a pre-computed grid (if the instrument/filter
combination is available in the grid) or at runtime (if the filter throughput
can be calculated but is not present in the grid). A Phoenix table
specifies each source with the following columns:

* ID (used to report the source in output tables)

* Dataset (all sources within a dataset must have the same age and metallicity)

* RA (decimal degrees)

* DEC (decimal degrees)

* Age (years)

* Metallicity ([Fe/H])

* Distance (kpc)

* Mass (solar masses)

* Binary status (whether or not the source is a binary companion of the previous source)

In order to identify a table as a Phoenix table, the following metadata must be
present:

* ``type=phoenix``

Phoenix tables will be converted to the Internal input format by the STIPS
Observation Module, and observed sources will be listed in an output table.

.. _realtime-phoenix-tables:

Realtime Phoenix Tables
-----------------------

This is a variant of Phoenix tables in which, even if a pre-calculated source
grid is present, it will never be used, and all stellar spectra will be
generated and observed at runtime. As a result, Realtime Phoenix tables are much
slower to generate. A Realtime Phoenix table specifies each source with the
following columns:

* ID

* RA (decimal degrees)

* DEC (decimal degrees)

* Teff

* Log(g)

* Metallicity ([Fe/H])

* Apparent Magnitude

In order to use a Realtime Phoenix table, the following metadata must be
present:

* ``type=phoenix_realtime``

* ``bandpass`` (this must be set to the bandpass in which the effective magnitude
  is measured)

Realtime Phoenix tables will be converted to the Internal input format by the
STIPS Observation Module, and observed sources will be listed in an output
table.

.. _pandeia-tables:

Pandeia Tables
--------------

Pandeia tables are usually used for internal testing of STIPS against Pandeia_
in order to ensure that STIPS results are sufficiently close to the results
produced by Pandeia. They are a variant of the Realtime Phoenix table, in which
the following columns are present:

* ID

* RA (decimal degrees)

* DEC (decimal degrees)

* Key

* Apparent Magnitude

Here the Key column replaces effective temperature, log(g), and metallicity,
and it is set to the key that Pandeia uses in producing its own Phoenix model
spectra. As such, it provides a more compact interface than a Realtime Phoenix
table at the cost of only allowing sources that are available as pre-defined
keys in Pandeia.

In order to use a Pandeia table, the following metadata must be present:

* ``type=pandeia``

* ``bandpass``

Bandpass is treated as it is in Realtime Phoenix tables. Table conversions are
done in exactly the same way as Realtime Phoenix tables.

.. _bc95-catalog:

BC95 Catalog
------------

A BC95 catalog is intended to include galaxies created from the `Bruzual and
Charlot Isochrone Synthesis Spectral Evolutionary Code (December 1995 version)
<https://stsynphot.readthedocs.io/en/latest/stsynphot/appendixa.html#bruzual-charlot-atlas>`_.
A BC95 catalog is an extended-source catalog, and specifies sources with the
following columns:

* ID

* RA (decimal degrees)

* DEC (decimal degrees)

* Redshift

* Model (one of 'a', 'b', 'c', 'd', or 'e', with the description of each model
  provided in the
  `BC95 README <https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/the-bruzual-charlot-atlas>`_
  file.

* Age (one of 10E5, 25E5, 50E5, 76E5, 10E6, 25E6, 50E6, 10E7, 50E7, 10E8, 50E8,
  10E9, years)

* Profile (one of 'expdisk' or 'devauc')

* Radius (arcseconds)

* Axial Ratio

* PA (degrees)

* Apparent Surface Brightness

In order to identify the catalog as a BC95 catalog, the following metadata
must be present:

* ``type=bc95``

* ``bandpass``

During the observation, the catalog will be converted into an Internal format,
with any necessary additional metadata added at this point. Galaxy spectra will
be generated from the atlas, and count rates derived through synphot observation
of the generated spectrum. An output catalog will be generated showing the
observed sources (along with their Sersic profile data).

.. _internal-catalog:

Internal Catalog
----------------

An Internal catalog is intended to include either point or extended sources,
but is limited to a single filter. It must contain the following columns:

* ID

* RA (decimal degrees)

* DEC (decimal degrees)

* FLUX (for point sources, count rate in the specified filter, counts/s. For
  sersic profiles, surface brightness inside Re in the specified filter,
  counts/s)

* TYPE (either 'point' or 'sersic')

* N (Sersic profile index if TYPE is 'sersic', otherwise ignored)

* Re (half-light radius in pixels if TYPE is 'sersic', otherwise ignored)

* Phi (angle of PA in degrees if TYPE is 'sersic', otherwise ignored)

* Ratio (axial ratio if TYPE is 'sersic', otherwise ignored)

* Notes (any notes that are needed. Not used directly, but any notes will be
  retained in the observed catalog produced during the observation.)

In order to identify the catalog as an Internal catalog, and in order to use
it for STIPS observations, the following columns must be present:

* ``type=internal``

* ``filter``

``filter`` is the filter to which the catalog has been calibrated. This
catalog type will not be converted during observation, but an observed source
catalog will be generated.

.. _mixed-catalog:

Mixed Catalog
-------------

A Mixed catalog is identical to an Internal catalog, except that it
contains one additional column:

* Units (one of 'p' for photons/s, 'e' for electrons/s, 'j' for Jansky, or 'c'
  for counts/s.)

In order to identify the catalog as a Mixed catalog, the following metadata
must be present:

* ``type=mixed``

* ``filter``

This catalog will have its flux values converted to counts/s, and will then be
treated as an Internal catalog.

.. _multifilter-catalog:

Multifilter Catalog
-------------------

A Multifilter catalog is identical to an Internal catalog, except that it
does not have a filter specified in its metadata and, instead of having a FLUX
column, it has one or more columns, each named after an available filter, that
provide the source count rate in that filter.

A Multifilter catalog must have the following metadata:

* ``type=multifilter``

The appropriate filter's count rate will be renamed as 'flux' as the catalog
is converted to Internal format.

.. _generic-catalog:

Generic Catalog
---------------

A Generic catalog is a point-source catalog with the following columns:

* RA (decimal degrees)

* DEC (decimal degrees)

* One column for each desired filter, showing the count rate in that filter.

* (Optional) an ID column for each source.

No specific metadata is required.

.. note::

	If a ``type`` metadata field is present in a Generic catalog, it must not
	have any of the above values. If it does, the catalog will be treated as
	whatever catalog type its type field indicates, and will probably fail to
	process.

.. _IPAC: https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html
.. _BinTable: https://docs.astropy.org/en/stable/io/fits/#working-with-table-data
.. _Pandeia: https://jwst.etc.stsci.edu
