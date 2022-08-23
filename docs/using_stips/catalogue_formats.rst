STIPS Catalog Formats
=====================
.. note::

    If you use the STIPS Scene Module to create your observation scenes, then
    created catalogues will already be in supported formats (specifically the
    Phoenix format for stars and the BC95 format for galaxies). Alternately, if
    you create new STIPS catalogues by editing existing catalogues, then, as long
    as you keep to the existing format, any created catalogues should work.

General Formatting
------------------

In general, STIPS catalogues may be in one of two formats, IPAC_ text or FITS
BinTable_. Both are accessed through the `astropy table
<https://docs.astropy.org/en/stable/table/>`_ API. In the case of FITS tables,
STIPS will internally divide tables into chunks of 100,000 entries, putting each
chunk into its own extension.

In IPAC_ tables, metadata is presented at the start of the table, in lines with
the format::

    \key=value

In FITS tables, metadata is included as header keywords.

STIPS Table Formats
-------------------

  \• Phoenix

  \• Realtime Phoenix

  \• Pandeia

  \• BC95

  \• Internal

  \• Mixed

  \• Multifilter

  \• Generic

  \• Output

Phoenix Tables
--------------

Phoenix tables are tables of point sources that can be generated as models using
the `Phoenix <http://phoenix.ens-lyon.fr>`_ star simulator. STIPS uses
`stsynphot <https://stsynphot.readthedocs.io/en/latest/>`_ to generate its
phoenix spectra, either using a pre-computed grid (if the instrument/filter
combination is available in the grid) or at runtime (if the filter throughput
can be calculated but is not present in the grid). A phoenix table
specifies each source with the following columns:

  \• ID (used to report the source in output tables)

  \• Dataset (all sources within a dataset must have the same age and metallicity)

  \• RA (decimal degrees)

  \• DEC (decimal degrees)

  \• Age (years)

  \• Metallicity ([Fe/H])

  \• Distance (kpc)

  \• Mass (solar masses)

  \• Binary status (whether or not the source is a binary companion of the previous
  source)

In order to identify a table as a Phoenix table, the following metadata must be
present:

  \• type=phoenix

Phoenix tables will be converted to the Internal input format by the STIPS
Observation Module, and observed sources will be listed in an output table.

Realtime Phoenix Tables
-----------------------

This is a variant of phoenix tables in which, even if a pre-calculated source
grid is present, it will never be used, and all stellar spectra will be
generated and observed at runtime. As a result, realtime phoenix tables are much
slower to generate. A realtime phoenix table specifies each source with the
following columns:

  \• ID

  \• RA (decimal degrees)

  \• DEC (decimal degrees)

  \• Teff

  \• Log(g)

  \• Metallicity ([Fe/H])

  \• Apparent Magnitude

In order to use a realtime phoenix table, the following metadata must be
present:

  \• ``type=phoenix_realtime``

  \• ``bandpass`` (this must be set to the bandpass in which the effective magnitude
  is measured)

Realtime phoenix tables will be converted to the Internal input format by the
STIPS Observation Module, and observed sources will be listed in an output
table.

Pandeia Tables
--------------

Pandeia tables are usually used for internal testing of STIPS against Pandeia_
in order to ensure that STIPS results are sufficiently close to the results
produced by Pandeia. They are a variant of the Realtime phoenix table, in which
the following columns are present:

  \• ID

  \• RA (decimal degrees)

  \• DEC (decimal degrees)

  \• Key

  \• Apparent Magnitude

Here the Key column replaces effective temperature, log(g), and metallicity,
and it is set to the key that pandeia uses in producing its own phoenix model
spectra. As such, it provides a more compact interface than a realtime phoenix
table at the cost of only allowing sources that are available as pre-defined
keys in pandeia.

In order to use a pandeia table, the following metadata must be present:

  \• type=pandeia

  \• bandpass

Bandpass is treated as it is in Realtime phoenix tables. Table conversions are
done in exactly the same way as Realtime phoenix tables.

BC95 Catalogue
--------------

A BC95 catalogue is intended to include galaxies created from the `Bruzual and
Charlot Isochrone Synthesis Spectral Evolutionary Code (December 1995 version)
<https://ssb.stsci.edu/pysynphot/docs/appendixa.html#pysynphot-appendixa-bc95">`_.
A BC95 catalogue is an extended-source catalogue, and specifies sources with the
following columns:

  \• ID

  \• RA (decimal degrees)

  \• DEC (decimal degrees)

  \• Redshift

  \• Model (one of 'a', 'b', 'c', 'd', or 'e', with the description of each model
  provided in the
  `BC95 README <https://www.stsci.edu/hst/observatory/crds/cdbs_bc95.html>`_
  file.

  \• Age (one of 10E5, 25E5, 50E5, 76E5, 10E6, 25E6, 50E6, 10E7, 50E7, 10E8, 50E8,
  10E9, years)

  \• Profile (one of 'expdisk' or 'devauc')

  \• Radius (arcseconds)

  \• Axial Ratio

  \• PA (degrees)

  \• Apparent Surface Brightness

In order to identify the catalogue as a bc95 catalogue, the following metadata
must be present:

  \• ``type=bc95``

  \• ``bandpass``

During the observation, the catalogue will be converted into an internal format,
with any necessary additional metadata added at this point. Galaxy spectra will
be generated from the atlas, and count rates derived through synphot observation
of the generated spectrum. An output catalogue will be generated showing the
observed sources (along with their sersic profile data).

Internal Catalogue
------------------

An Internal catalogue is intended to include either point or extended sources,
but is limited to a single filter. It must contain the following columns:

  \• ID

  \• RA (decimal degrees)

  \• DEC (decimal degrees)

  \• FLUX (for point sources, count rate in the specified filter, counts/s. For
  sersic profiles, surface brightness inside Re in the specified filter,
  counts/s)

  \• TYPE (either 'point' or 'sersic')

  \• N (Sersic profile index if TYPE is 'sersic', otherwise ignored)

  \• Re (half-light radius in pixels if TYPE is 'sersic', otherwise ignored)

  \• Phi (angle of PA in degrees if TYPE is 'sersic', otherwise ignored)

  \• Ratio (axial ratio if TYPE is 'sersic', otherwise ignored)

  \• Notes (any notes that are needed. Not used directly, but any notes will be
  retained in the observed catalogue produced during the observation.)

In order to identify the catalogue as an internal catalogue, and in order to use
it for STIPS observations, the following columns must be present:

  \• ``type=internal``

  \• ``filter``

``filter`` is the filter the catalogue has been calibrated to. This catalogue type
will not be converted during observation, but an observed source catalogue will
be generated.

Mixed Catalogue
---------------

A Mixed catalogue is identical to an internal catalogue, except that it
contains one additional column:

  \• Units (one of 'p' for photons/s, 'e' for electrons/s, 'j' for Jansky, or 'c'
  for counts/s.)

In order to identify the catalogue as a mixed catalogue, the following metadata
must be present:

  \• ``type=mixed``

  \• ``filter``

This catalogue will have its flux values converted to counts/s, and will then be
treated as an internal catalogue.

Multifilter Catalogue
---------------------

A Multifilter catalogue is identical to an internal catalogue, except that it
does not have a filter specified in its metadata and, instead of having a FLUX
column, it has one or more columns, each named after an available filter, that
provide the source count rate in that filter.

A Multifilter catalogue must have the following metadata:

  \• ``type=multifilter``

The appropriate filter's count rate will be renamed as 'flux' as the catalogue
is converted to internal format.

Generic Catalogue
-----------------

A Generic catalogue is a point-source catalogue with the following columns:

  \• RA (decimal degrees)

  \• DEC (decimal degrees)

  \• One column for each desired filter, showing the count rate in that filter.

  \• (Optional) an ID column for each source.

No specific metadata is required.

.. note::

	If a 'type' metadata field is present in a generic catalogue, it must not
	have any of the above values. If it does, the catalogue will be treated as
	whatever catalogue type its type field indicates, and will probably fail to
	process.

.. _IPAC: https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html
.. _BinTable: https://docs.astropy.org/en/stable/io/fits/#working-with-table-data
.. _Pandeia: https://jwst.etc.stsci.edu
