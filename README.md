Pipeline to reduce NACO L' (AGPM and satPSF) and M' data

Prerquisits
===========

-target_properties.pro
  -first create a file target_ids.txt containing one column with the file ID. Has to be identical to the target ID in the fits header and should not contain spaces!
  -after changing the path parameters (line 18) run target_properties. it creates an IDL .sav file for each target containing stellar properties as derived from simbad

-file structure
  -each NACO observation has to be in a unique directory, <target ID>/RAW/
  -the target ID has to be identical to the ID provided in the fits header and has to be identical to the .sav file! If there are multiple epochs it is OK to use, e.g. <target ID>_1stEpoch/RAW or similar.
  -If no AGPM was used the directory has to be named <target ID>_noAGPM/RAW/
  -there has to be a text file "datapaths.txt" containing one row with the full path to the individiual RAW directory, e.g. /path/to/the/directory/HD123/ (without including "RAW")
  
  
-the RAW directory of each star contains all raw data, i.e. science and calibration files:
  -all science frames (ADI sequence, unsaturated flux)
  -darks for science and unsaturated flux measurements
  -flat fields (5 images)
  -there can be only files used during the reduction process, i.e. no multiple flat field sequences, no multiple dark images from different days, no science cubes with different DITs, no unsaturated flux files with different DITs, ...only cleaned raw data
  

Reduction
=========

-adjust path in batch_NACO_reduce.pro for datapaths.txt (line 107)
-adjust path in batch_NACO_reduce.pro for the location of the .sav files produced by target_properties.pro (line 124)
-run batch_NACO_reduce.pro and follow instructions
-final products will be img_<filter>_dc.fits, vec_<filter>_paral.fits, PSF_<filter>.fits, where <filter> is "Lp" or "Mp" stored in <target ID>/Reduced/
