
import numpy as np
import matplotlib.pyplot as plt

# 1D spherical symmetry diffusion: dC/dt = D*(d2C/dr2 + 2/r*dC/dr) - kdeg*C + S(r,t)
# Explicit scheme; ensure stability: dt < dr^2/(2*D)
D = 5e-6   # mm^2/s  (adjust for medium)
kdeg = 0.0 # 1/s
R = 0.5    # mm radius of spheroid/domain
Nr = 200
dr = R/Nr
T_s = 6*3600  # total time (s)
dt = 0.2*dr**2/D  # conservative stability step
Nt = int(T_s/dt)

# Source: bacteria density localized near core or shell
def rho_bac(r):
    # example: uniform density inside radius 0.4R
    return 1.0 if r <= 0.4*R else 0.2

q_glu = 1e-7  # mM/s per density unit (placeholder)

r = np.linspace(0, R, Nr+1)
C = np.zeros_like(r)
C_new = C.copy()

for n in range(Nt):
    for i in range(1, Nr):
        rp = r[i]+1e-9
        d2 = (C[i+1] - 2*C[i] + C[i-1]) / dr**2
        d1 = (C[i+1] - C[i-1]) / (2*dr)
        lap = d2 + 2.0/rp * d1
        S = q_glu * rho_bac(r[i])
        C_new[i] = C[i] + dt*(D*lap - kdeg*C[i] + S)
    # Neumann BC at r=0 (symmetry): dC/dr=0 -> C[-1] mirrored
    C_new[0] = C_new[1]
    # Dirichlet or Robin at boundary r=R; here Dirichlet C=0 (sink) as example
    C_new[-1] = 0.0
    C, C_new = C_new, C

# Save profile
np.savetxt("spheroid_glu_profile.csv", np.c_[r, C], delimiter=",", header="radius_mm,glu_mM", comments="")
plt.figure(); plt.plot(r, C); plt.xlabel("radius (mm)"); plt.ylabel("Glu (mM)"); plt.title("End-time Glu profile"); plt.savefig("spheroid_glu_profile.png")
print("Saved spheroid_glu_profile.csv and spheroid_glu_profile.png")
