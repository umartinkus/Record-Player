"""Full-chain system model: cartridge mechanics -> cartridge electrical
network -> RIAA amplifier.

Python port of system_model.m.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from cartridge_tf_elec import cartridge_tf_elec
from cartridge_tf_mech import cartridge_tf_mech
from riaa_amp_tf import riaa_amp_tf

# L1_mH = 0.80e-3; C2_uF = 25e-6; C1_uF = 0.0506e-6; R_ohm = 40
# L1_mH = 0.8e-3; C1_uF = 0.0506e-6; C2_uF = 17e-6; R_ohm = 40
# Hvel, Hforce = cartridge_tf_mech(L1_mH, C1_uF, C2_uF, R_ohm)
#
# print(f'Series (force notch):  {1/(2*np.pi*np.sqrt(L1_mH*C2_uF)):.0f} Hz')
# print(f'Parallel (force peak): {1/(2*np.pi*np.sqrt(L1_mH*C1_uF)):.0f} Hz')
#
# (bode of Hvel and Hforce from 2*pi*20 to 2*pi*50e3 rad/s)


def series(*systems):
    """Cascade TransferFunctions (MATLAB's sys1*sys2*...)."""
    num, den = np.array([1.0]), np.array([1.0])
    for sys in systems:
        num = np.polymul(num, sys.num)
        den = np.polymul(den, sys.den)
    return signal.TransferFunction(num, den)


Hvel, _ = cartridge_tf_mech()  # MATLAB's single-output call returns Hvel
Hel, _ = cartridge_tf_elec()
sys = series(Hvel, Hel, riaa_amp_tf())

w = 2 * np.pi * np.logspace(0, 5, 1000)  # 1 Hz to 100 kHz
w, mag_dB, phase_deg = signal.bode(sys, w)

fig, (ax_mag, ax_ph) = plt.subplots(2, 1, sharex=True)
ax_mag.semilogx(w, mag_dB)
ax_mag.set_ylabel('Magnitude (dB)')
ax_mag.set_title('Bode Diagram')
ax_mag.grid(True, which='both')
ax_ph.semilogx(w, phase_deg)
ax_ph.set_ylabel('Phase (deg)')
ax_ph.set_xlabel('Frequency (rad/s)')
ax_ph.grid(True, which='both')

plt.show()
