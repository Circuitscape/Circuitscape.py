16/12/11
========

* Check for conditioning of the input by looking for ratios of largest and
  smallest values not being larger than a threshold such as 10000.

* Investigate the use of MKL.

* Investigate the use of OpenMP to speed up matvecs in scipy.


12/7/11
=======

* Add network mode checkbox, and add pairwise network to list of
  modeling modes.

* Add advanced network mode?  Would need to calculate cumulative
  voltages, would need to do components for graph, sources, and
  grounds.  Could use focal node component code for latter.

* Average current across core area? Probably not.

* Add easier input for exluded pairs/variable soruces- list instead of
   matrix format.

* Add include_exclude pairs to network (low priority)

* Remove log transform option? 


3/14/09
=======

* Reorganize cs_compute. Break into smaller files - single_ground,
  multiple_ground, solver, file i/o, logging etc.

* The writing of current and voltage maps should be done inside the
  create_voltage_map and create_current_map routines, rather than in
  the single_ground and multiple_ground modules.

* Decouple the GUI and cs_compute further. Allow for possibility of
  GUI and compute to run on different computers in the future.

* Plotting maps in the GUI.

* Rework the current map creation code to create all current maps in
  the end. Make it faster also.

* Improved logging information on the console. 
   
   a. Report a banner with version number, OS, python version and
      versions of numpy, scipy, pyamg, gdal and other libraries.
   
   b. Report problem stats - number of components, size of the graph,
      total numbers of solves.

* Use GDAL for file I/O. GDAL binaries are available on all
  platforms. Bundling may need some work.

* Integrate with ArcGIS and any other third party software of interest.

* Try to avoid doing A+A' in Laplacian. Perhaps the graph builder can be reformulated.

* Check 64-bit capability in cs_verify.

* Provide time and space requirements for scalable synthetic problems.

* Try umfpack for smaller problems.
