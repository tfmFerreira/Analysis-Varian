#converts VNMR data from 1D experiments into 3 columns (time, FID real, FID imag)
#Here loops over different folders that have identical names except for a number

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import os

rootdir = os.getcwd()

# The list (vector) of folder numbers you want to loop over
folder_numbers = [40754, 40754, 43000, 43916, 44101, 44131, 44500, 46775, 50297]

for num in folder_numbers:
	folder_name = f"20251031_13C_rinept_37C_DECpwr{num}.fid"
	dirname = os.path.join(rootdir, folder_name)
	dic,data = ng.varian.read(dirname,"fid","procpar")
	A = ng.varian.create_data(data)
	Bf1=dic.get("procpar").get("H1reffrq").get("values")
	NP = float(dic.get("procpar").get("np").get("values")[0])
	SW = float(dic.get("procpar").get("sw").get("values")[0])
	DW = 1/(SW)
	ACQtime = np.arange(NP/2) * DW
	out_file = os.path.join(dirname, f"output.txt")
	np.savetxt(out_file, np.column_stack((ACQtime, A.real, A.imag)), fmt="%.6f")
	Bf1 = np.array(Bf1, dtype=float)   # convert everything to floats
	out_file = os.path.join(dirname, f"Bf1.txt")
	np.savetxt(out_file, Bf1, fmt="%.6f")

