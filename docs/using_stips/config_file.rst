The STIPS Configuration File
============================
.. note::

    Although the configuration file is the recommended general way of setting STIPS options,
    it is also possible to set many of these options with dictionary keys or by specifying
    keyword arguments to classes or functions at runtime.

Configuration File Format
-------------------------

STIPS configuration files are formatted as `YAML files <https://yaml.org>`_. The
configuration files included in the internal STIPS data directory (or the
downloadable ``stips_data`` directory) include comments describing how the keywords
are used, and what values they can have. This commentary is optional.

Configuration File Strategy
---------------------------

When the ``stips.utilities.SelectParameter()`` function is used to look for a
configuration parameter, STIPS follows the following strategy to find that
parameter (and will take whichever value it finds first). As such, if the user
wishes to override only some STIPS configuration values, they can put a
configuration file with only those values specified into the directory where
STIPS is being run. Any other configuration will fall back to other files.

The STIPS configuration hierarchy is:

#. Keyword arguments provided when creating a STIPS class.
#. Any configuration file (or directory containing a file named
   ``stips_config.yaml``) provided directly to
   ``stips.utilities.SelectParameter()``.
#. A file named ``stips_config.yaml`` in the same directory as STIPS is being run.
#. A file (or a directory containing a file named ``stips_config.yaml``)
   named in the environment variable ``stips_config``.
#. A file named ``stips_config.yaml`` located in the ``stips_data`` directory.
#. The ``stips_config.yaml`` file in the internal STIPS data directory.

.. note::

	The internal ``stips_config.yaml`` file should not be edited, as it contains
	default values that should not be changed directly. The appropriate way for
	a user to change the STIPS configuration is to use any of methods 1-5 above.


General Configuration Keywords
------------------------------

input_location (default *$CWD*)
	If classes are given an input file name which is not a fully-qualified path,
	they will look for that file in this directory. ``$CWD`` is a special value
	that is replaced at runtime with the directory from which STIPS is being
	run. If used as a keyword argument, ``cat_path`` can be used instead of
	``input_location`` as a variable name for historical reasons.

output_location (default *$CWD$*)
	If classes produce any output files, they will be placed in this directory.
	``$CWD`` is a special value that is replaced at runtime with the directory
	from which STIPS is being run. If used as a keyword argument, ``out_path``
	can be used instead of ``output_location`` for historical reasons.

catalogue_type (default *fits*)
	This is the format in which catalogue files will be written, and from which
	they will be read by default. Currently the only supported values for this
	keyword are ``ascii.ipac`` and ``fits``. If used as a keyword argument,
	``cat_type`` can be used instead of ``catalogue_type`` for historical reasons.

log_level (default *INFO*)
	This is the default log level for the internal STIPS logger. Any value that
	can be obtained from the logging module using ``getattr(logging, VALUE)``
	and then passed to ``logging.logger.setLevel()`` can be used here.

random_seed (default *1234*)
	The seed value to be passed to ``numpy.random.RandomState(seed=seed)``. A
	value of "-1" means that seed will be treated as None (and so will use
	/dev/urandom or the system clock). If used as a keyword argument,
	``seed`` can be used instead of ``random_seed`` for historical reasons.


Observation Keywords
--------------------

observation_default_background (default *0.0*)
	The default sky background in counts/s/detector pixel. Currently this keyword can be set to:

	* any integer or floating point value, in which case that value will be
	  used directly.

	* any of the string values 'none', 'low', 'avg', or 'high'. In this
	  case, 'none' is always treated as zero, and for any other keyword if the
	  value is defined for the instrument/detector selected, that value will
	  be used. If no such value can be found, the background will be set to 0.

	If used as a keyword argument, ``background`` can be used instead of
	``observation_default_background`` for historical reasons.

observation_jbt_location (default *$WEB*)
	If JBT is being used to determine the background, this tells STIPS where the
	JBT data is located. ``$WEB`` indicates that the value should be fetched
	from online, ``$DATA`` indicates that the value should be taken from a
	directory named ``background`` in the ``stips_data`` directory. Otherwise,
	the value should be the path to a directory containing a local cache of the
	data. If used as a keyword argument, ``background_location`` or
	``jbt_location`` can be used instead of ``observation_jbt_location`` for
	historical reasons.

observation_distortion_enable (default *false*)
	Whether co-ordinate distortion information should be included in the
	observation. Note that this is not yet available for Roman. If used as a
	keyword argument, ``distortion`` can be used instead of
	``observation_distortion_enable`` for historical reasons.


PSF Convolution Configuration
-----------------------------

psf_grid_default_size (default *1*)
	What size PSF grid should be created. Note that this value is expressed as
	a side length, so if psf_grid_default_size is set to n, WebbPSF will create
	a total of n^2 PSF images. If used as a keyword argument, ``psf_grid_size``
	can be used instead of ``psf_grid_default_size`` for historical reasons.

psf_cache_enable (default *true*)
	Whether PSF grids created by webbpsf should be cached after creation for
	potential re-use.

psf_cache_location (default *$DATA*)
	Where PSF grids should be cached if caching in enabled. The special value
	``$DATA`` indicates the stips_data directory.

psf_cache_directory (default *psf_cache*)
	The name of the directory inside ``psf_cache_location`` where PSF grids
	should be cached (again, if caching is enabled).

psf_convolution_max_size (default *8192*)
	The maximum data array size to create for convolutions. Note that this value
	should be a power of two. If it isn't, the largest power of two that is less
	than ``psf_convolution_max_size`` will be used. If used as a keyword
	argument, ``convolve_size`` can be used instead of
	``psf_convolution_max_size`` for historical reasons.


Error Residual Configuration
----------------------------

residual_convolve_psf (default *true*)
	Whether PSF convolution should be performed when adding error. If used as a
	keyword argument, ``convolve`` may be used instead of ``residual_convolve_psf``
	for historical reasons.

residual_poisson (default *true*)
	Whether Poisson noise should be added when adding error.

residual_readnoise (default *true*)
	Whether Readnoise should be added when adding error.

residual_flat (default *true*)
	Whether a flatfield removal residual should be added when adding error.

residual_dark (default *true*)
	Whether a dark current removal residual should be added when adding error.

residual_cosmic (default *true*)
	Whether cosmic ray removal residuals should be added when adding error.
