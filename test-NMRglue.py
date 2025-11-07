import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import os

dirname = os.getcwd()

print(dirname)

dic,data = ng.varian.read(dirname,"fid","procpar")

A = ng.varian.create_data(data)

Bf1=dic.get("procpar").get("H1reffrq").get("values")
NP = float(dic.get("procpar").get("np").get("values")[0])
SW = float(dic.get("procpar").get("sw").get("values")[0])
DW = 1/(SW)
ACQtime = np.arange(NP/2) * DW
np.savetxt("output.txt", np.column_stack((ACQtime, A.real, A.imag)), fmt="%.6f")
Bf1 = np.array(Bf1, dtype=float)   # convert everything to floats
np.savetxt("Bf1.txt", Bf1, fmt="%.6f")

