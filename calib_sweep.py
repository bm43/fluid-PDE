# calibrate_and_sweep.py
# Calibrate to Goldmann and sweep permeability/EVP/inflow following ocular practice [web:225][web:1]
import numpy as np

# Given TM thickness L, area A, viscosity mu, Darcy permeability k:
# Facility C ≈ (k * A) / (mu * L) in SI, linking Darcy resistance to Goldmann outflow facility [web:225]
mu = 0.0007
A_tm = 2e-4 * 2.5e-3  # height * circumferential wedge width proxy [web:235]
L_tm = 2e-4

def facility_from_perm(k):
    return (k * A_tm) / (mu * L_tm)

def goldmann_iop(EVP, Q, C):
    # IOP ≈ EVP + Q/C (ignoring unconventional outflow for baseline) [web:225]
    return EVP + Q / C

# Sweep examples
EVP_vals = np.array([1200, 1500, 1800])  # Pa (~9, 11, 13 mmHg) [web:225]
Q = 2.5e-9/60.0                           # m^3/s inflow [web:6]
perms = np.logspace(-15, -13, 5)          # TM permeability sweep [web:1]

rows = []
for evp in EVP_vals:
    for k in perms:
        C = facility_from_perm(k)
        iop = goldmann_iop(evp, Q, C)
        rows.append((evp, k, C, iop))

# Save calibration table
import csv
with open("goldmann_sweep.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["EVP_Pa", "k_tm_m2", "facility_m3sPa", "IOP_Pa"])
    w.writerows(rows)
