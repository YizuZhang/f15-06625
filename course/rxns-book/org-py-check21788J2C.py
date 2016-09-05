import numpy as np
from scipy.integrate import odeint

vo = 476.0    # ft^3 / hr
Po = 35.0     # atm
T = 1500.0    # Rankine
R = 0.73      # in appropriate units
CTo = Po / R / T
Fto = CTo * vo

# initial molar flows
Fmo = Fto / 3.0
Fho = Fto * 2.0 / 3.0
Fxo = 0.0
Fmeo = 0.0
Ftolo = 0.0

def dFdV(F, t):
    'PFR moe balances'
    Ft = F.sum()

    v = vo * Ft / Fto
    C = F / v
    CM, CH, CX, CMe, CT = C

    # rate laws
    k1 = 55.20
    k2 = 30.20
    r1m = -k1 * CM * CH**0.5
    r2t = k2 * CX * CH**0.5

    # net rates
    rM = r1m
    rH = r1m - r2t
    rX = -r1m - r2t
    rMe = -r1m + r2t
    rT = r2t

    dFMdV = rM
    dFHdV = rH
    dFXdV = rX
    dFMedV = rMe
    dFTdV = rT

    return [ dFMdV, dFHdV, dFXdV, dFMedV, dFTdV ]

Finit = [Fmo, Fho, Fxo, Fmeo, Ftolo]
Vspan = np.linspace(0.0, 238.0)

sol = odeint(dFdV, Finit, Vspan)

Ft = sol.sum(axis=1)  # sum each row
v = vo * Ft / Fto

FM  = sol[:,0]
FH  = sol[:,1]
FX  = sol[:,2]
FMe = sol[:,3]
FT  = sol[:,4]
F1,F2,F3,F4,F5 = [sol[:,i] for i in range(5)]
F1,F2,F3,F4,F5 = sol.T
F1,F2,F3,F4,F5 =  np.transpose(sol)

tau = Vspan / vo

import matplotlib.pyplot as plt
plt.plot(tau, FM / v, label='$C_M$')
plt.plot(tau, FH / v, label='$C_H$')
plt.plot(tau, FX / v, label='$C_X$')

plt.legend(loc='best')
plt.xlabel('$\\tau$ (hr)')
plt.ylabel('Concentration (lbmol/ft$^3$)')
plt.savefig('images/multiple-reactions-pfr.png')
plt.show()
