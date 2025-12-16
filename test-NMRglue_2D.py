#converts VNMR data from 2D experiments into x columns (time, FID real, FID imag)
#Here loops over different folders that have identical names except for a number


import os
import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np

rootdir = os.getcwd()

# The list (vector) of folder numbers you want to loop over
folder_numbers = np.arange(1, 10)

for num in folder_numbers:
    folder_name = f"20251212_9_13C_B1inhomog_{num}.fid"
    dirname = os.path.join(rootdir, folder_name)
    dic,data = ng.varian.read(dirname,"fid","procpar")
    A = ng.varian.create_data(data)
    Bf1=dic.get("procpar").get("H1reffrq").get("values")
    NP = float(dic.get("procpar").get("np").get("values")[0])
    SW = float(dic.get("procpar").get("sw").get("values")[0])
    DW = 1/(SW)
    ACQtime = np.arange(NP/2) * DW
    out_file_real = os.path.join(dirname, f"output_real.txt")
    out_file_imag = os.path.join(dirname, f"output_imag.txt")
    data_real = np.column_stack([ACQtime] + [A[i, :].real for i in range(A.shape[0])])
    data_imag = np.column_stack([ACQtime] + [A[i, :].imag for i in range(A.shape[0])])
    np.savetxt(out_file_real, data_real, fmt="%.6f")
    np.savetxt(out_file_imag, data_imag, fmt="%.6f")
    Bf1 = np.array(Bf1, dtype=float)   # convert everything to floats
    out_file = os.path.join(dirname, f"Bf1.txt")
    np.savetxt(out_file, Bf1, fmt="%.6f")
    increment=dic.get("procpar").get("pwXstep").get("values")
    increment = np.array(increment, dtype=float)   # convert everything to floats
    out_file = os.path.join(dirname, f"increment.txt")
    np.savetxt(out_file,increment, fmt="%.6f")
    coarsepwr=float(dic.get("procpar").get("tpwr").get("values")[0])
    finepwr=float(dic.get("procpar").get("aX90").get("values")[0])
    out_file = os.path.join(dirname, f"power.txt")
    with open(out_file, "w") as f:
        f.write(f"{coarsepwr} {finepwr}\n")
