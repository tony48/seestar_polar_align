from cx_Freeze import setup, Executable
import sys

# Dependencies are automatically detected, but it might need
# fine tuning.

buildOptions = dict(include_files = [('dwarf_ble_connect/','./dwarf_ble_connect'),('Install/','.')]) 
#folder,relative path. Use tuple like in the single file to set a absolute path.

 
base = 'Win32GUI' if sys.platform=='win32' else None
setup(
    name = "PhotoPolarAlign",
    version = "1.0",
    description = "Photo Polar Align with the Dwarf 2",
    options = dict(build_exe = buildOptions),
    executables = [Executable("PPA.py")]
)