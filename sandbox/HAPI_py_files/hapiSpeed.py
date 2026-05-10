import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hapi import *
import time
import numpy as np
fetch('CO2',2,1,6000,6400)

# exit()
# start = time.time()
# nu_l, cs_co2_Voigt290 = absorptionCoefficient_Voigt(SourceTables='CO2',WavenumberStep=0.001, WavenumberRange=[6000,6400],Environment={'p':1,'T':293})

# doppler_times = []
# lorentz_times = []
# voigt_times = []
#
# for i in range(100):
#     start = time.time()
#     absorptionCoefficient_Voigt(SourceTables='CO2',WavenumberStep=0.001, WavenumberRange=[6000,6400],Environment={'p':1,'T':293})
#     voigt_times.append(time.time() - start)
#
#     print(voigt_times[-1])
#
# print("average voigt: %s" % (np.mean(np.array(voigt_times))))
# print("std dev voigt: %s" % (np.std(np.array(voigt_times))))

start = time.time()
nu, voigt_test_python = absorptionCoefficient_Voigt(SourceTables='CO2',WavenumberStep=0.001, WavenumberRange=[6000,6400],Environment={'p':1.0,'T':296.0}, WavenumberWing=40)
print(time.time() - start)

voigt_test_python = np.array(voigt_test_python)
voigt_test_julia = np.loadtxt("voigt_test_julia.csv")
print("Elapsed: %s s" % (time.time() - start))

# print(voigt_test_python.size)
# print(voigt_test_julia.size)
# print(np.arange(6000,6400.001,0.001).size)

# plt.plot(np.arange(6000,6400.0001,0.001), voigt_test_python, label="python")

plt.plot(np.linspace(6000,6400,400001), voigt_test_python, label="python")
plt.plot(np.linspace(6000,6400,400001), voigt_test_julia, label="julia")
plt.title("voigt Lineshape")
plt.xlabel("nu")
plt.legend()
plt.xlim(6300, 6350)
plt.savefig("voigt_pyjl.png")

plt.clf()

plt.plot(np.linspace(6000,6400,400001), np.divide((voigt_test_python-voigt_test_julia), voigt_test_python))
plt.title("Relative Difference: (voigt_python - voigt_julia) / voigt_python")
plt.xlabel("nu")
plt.xlim(6300, 6350)
plt.savefig("voigt_pyjl_relative.png")

plt.clf()

plt.plot(np.linspace(6000,6400,400001), voigt_test_python-voigt_test_julia)
plt.title("Absolute Difference: voigt_python - voigt_julia")
plt.xlabel("nu")
plt.xlim(6300, 6350)
plt.savefig("voigt_pyjl_absolute.png")
