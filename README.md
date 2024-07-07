# PhotoPolarAlign
A python utility to help align equatorial telescopes by imaging the Celestial Pole region.

Inspired by Dave Rowe http://www.considine.net/aplanatic/align.htm

Working with Python3

python3 -m pip install -r requirements

-----------------------------------------
New Version for Dwarf 2

Installation

 1. First Install libraries with : 
  
      python -m pip install -r requirements

 2. Then you need the dwarf_python_api found also on my github

      python -m pip install -r requirements-local.txt --target .

    This library must be installed locally in the root path of this project with the parameter --target .

 4. To generate the exe file from python, I'm using Cx Freeze, install it with:

       python -m pip install cx_Freeze
     
    Some extra files are need, to generate the exe: image, config files and some files from astropy
    that are utf8 but Cx Freeze copy them in ASCII so need this
    Unzip the Install Extra Files for Setup.zip, this create a Install directory.

    Now build the exe with:

      python setup.py build

    The final executable with all the files will be in the subdir of build directory.

 3. Configure the config.ini with your specifics data:

      [CONFIG]
      longitude = 
      latitude = 
      timezone = 
      local_photo_directory = .\Images
      ble_psd = DWARF_12345678
      ble_sta_ssid = 
      ble_sta_pwd = 

   The timezone value is like this

       timezone = Europe/Paris

   The last three corresponds to your home wifi values

 4. Finally you can start it with:

      python .\PPA.py
 
 5. To resolve your photos live you will need an api key on Nova Astrometry, setup in the settings Menu

      https://nova.astrometry.net/api_help