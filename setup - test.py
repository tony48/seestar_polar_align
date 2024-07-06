from cx_Freeze import setup, Executable
import os, sys
import glob


env_root = os.path.abspath(os.path.join(os.path.split(sys.executable)[0], ".."))
astropyPath = os.path.join(env_root, "Lib", "site-packages", "astropy", "*")

print(astropyPath)
fileList = glob.glob(astropyPath)
print(fileList)

moduleList = []

for mod in fileList:
    if os.path.isdir(mod):  # Check if the path is a directory
        modules = os.path.basename(mod)
        print(modules)
        moduleList.append("astropy." + modules)

print(moduleList)

# Dependencies are automatically detected, but it might need
# fine tuning.
build_options = {'packages': moduleList, 'excludes': [], 'zip_include_packages':''}
 
base = 'Win32GUI' if sys.platform=='win32' else None
setup(
    name = "PhotoPolarAlign",
    version = "0.1",
    description = "Photo Polar Align with the Dwarf 2",
    executables = [Executable("PPA.py")]
)