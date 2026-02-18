from hapi import *
import time
import numpy as np

names = ["H2O", "CO2", "O3", "N2O", "CO"]#, "CH4", "O2", "NO", "SO2", "NO2"]
molecules = [1, 2, 3, 4, 5]#, 6, 7, 8, 9, 10]
isotopes = [1, 1, 1, 1, 1]#, 2, 1, 2, 1, 1]

temperatures = [100, 175, 250, 325, 400]
pressures = [250, 500, 750, 1000, 1250]

fetch('CO2',2,1,6000,6400)

i = 1; 

for temp in temperatures:
    for pres in pressures:
        print i
        start = time.time()
        nu, voigt_test = absorptionCoefficient_Voigt(SourceTables='CO2',WavenumberStep=0.01, WavenumberRange=[6000,6400],Environment={'p':pres/1013.25,'T':temp}, WavenumberWing=40)
        print(time.time() - start)
        voigt_test = np.array(voigt_test)
        np.savetxt("Voigt_CO2_T" + str(temp) + "_P" + str(pres).replace('.', '') + ".csv", voigt_test, delimiter=',')
        i = i + 1

# P = 1000
# T = 250
# for i in range(len(names)):
#     print names[i]
#     fetch(names[i],molecules[i],isotopes[i],6000,6400)
#     start = time.time()
#     nu, voigt_test = absorptionCoefficient_Voigt(SourceTables=names[i],WavenumberStep=0.01, WavenumberRange=[6000,6400],Environment={'p':P/1013.25,'T':T}, WavenumberWing=40)
#     print(time.time() - start)
#     voigt_test = np.array(voigt_test)
#     np.savetxt("Voigt_" + names[i] + "_T" + str(T) + "_P" + str(P).replace('.', '') + ".csv", voigt_test, delimiter=',')