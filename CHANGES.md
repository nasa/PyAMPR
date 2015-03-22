Change Log
----------
v1.4 major changes (03/22/2015):
1. Added support for Level 2B netCDF files. Starting with IPHEX,
   processed AMPR instrument files will be provided in a netCDF-4 format. 
2. Renamed read_ampr_tb_leve1b to read_ampr_tb_level2b. Technically, the published
   AMPR data are Level 2B, and Level 1B data consist of pre-QC’d lower-level data. 
3. Swapped to standard AMPR convention for AmprTb.swath_angle and AmprTb.swath_left.
   Now go from -44.1 to 44.1 deg from left to right. Note AMPR scans right to left 
   but by convention reports its data left to right.

v1.3.2 major changes (09/24/2014):
1. Deleted write_ampr_kml method from AmprTb and all associated helper functions in 
   google_earth_tools.
2. Removed private but dead code.

v1.3.1 major changes:
1. Global constant variable renaming to conform better with PEP8 standards.
2. Converted to installable module (complete with setup.py script).
3. Added support for reading gzipped AMPR TB files without first decompressing

v1.3 major changes:
1. Significant refactoring of major methods to make use of new internal 
   methods. This reduced code duplication and makes the major methods
   more readable. It also reduced the number of local variables in the
   major methods.
2. Significant variable renaming to improve clarity. Most notably,
   Ampr_Tb class is now AmprTb to conform to PEP8 standard.
3. Put hard coding of constants at top of program to make it obvious.
4. Added the ability to suppress plotting of swaths during aircraft maneuvers
   in plot_ampr_track() and write_ampr_kmz()
5. Added timerange keyword to main plotting methods to allow filtering by
   time instead of scan number.
6. Added flag to enable returning of figure, Basemap, etc. objects from
   plot_ampr_track(). This empowers the creation of highly customized plots
   while still using plot_ampr_track() as a baseline for the AMPR component.
7. Added timespan stamp to write_ampr_kmz() and google_earth_tools.py. Now 
   multiple KMZ files from different times can be viewed in sequence in 
   Google Earth.

v1.2 major changes:
1. Further adjustments to the plotting routines to handle grossly 
   bad geolocations. Now, the presence of bad_data, 0s, -1s, or wildly
   varying Latitude/Longitude values within an individual scan gets the scan
   not considered for plotting purposes. There is a flag, equator, that can
   be set to True if the aircraft is flying near the Equator or Prime Meridian
   and 0s and -1s are normally expected.
2. Added internal methods _get_scan_indices() and _filter_bad_geolocations()
   to handle repeated tasks in the plotting routines.
3. Changed calc_polarization() to no longer consider aircraft roll angle.
   Now much quicker due to matrix multiplication. Flags are available 
   to switch between simple substitution or constrained linear inversion, 
   and to switch between forcing nadir matching in channels A & B or not.
   Can also now just call deconvolution for individual freqs if desired.
4. Changed default colormap to Brent Robert's amprTB_cmap, available in
   udf_cmap.py. If not available will just use cm.GMT_wysiwyg, available 
   from Basemap (a required dependency).
5. Adjusted the read routine to set infinites/NaNs in Level 1B TBs and 
   Latitude/Longitude to bad_data. Was blowing up Python during the
   plotting routines otherwise.

v1.1 major changes:
1. Added support for all projects on the GHRC, plus IPHEX.
2. Adjustments to the plotting routines to handle grossly bad geolocations
3. Adjusted plot_ampr_track() to put the colorbar on its own separate axis,
   so it is not varying in size with the plot.
4. Strip charts from plot_ampr_channels() now can show deconvolved H & V 
   via a keyword argument.
5. ER2 now called Aircraft to keep things generalized.
