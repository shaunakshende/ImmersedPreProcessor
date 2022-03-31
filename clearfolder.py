import glob, os
cwd=os.getcwd()
os.chdir(cwd)
for file in glob.glob("*.dat"):
    if file.startswith("XVOL") | file.startswith("input_coor") | file.startswith("Geometry") | file.startswith("foreground"):
        os.remove(file)
for file in glob.glob("*.txt"):
    os.remove(file)
