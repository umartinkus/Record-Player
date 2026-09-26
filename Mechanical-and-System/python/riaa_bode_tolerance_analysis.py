"""RIAA_BODE_TOLERANCE_ANALYSIS

Python port of riaa_bode_tolerance_analysis.m.

Bode plot of the RIAA amplifier (see riaa_amp_tf.py) with an
uncertainty band computed by first-order error propagation through
the transfer function, given assumed per-component tolerances.
Also overlays the ideal IEC-amended RIAA target curve for
reference, so you can see both "how far is the nominal design from
ideal" and "how much does part tolerance spread that."

METHOD (error propagation, not Monte Carlo):
At each frequency, the magnitude response in dB is a function of
six component values, dB = f(R0,C0,R1,C1,R2,C2). For each
component we estimate the local sensitivity dDB/dp_i numerically
(central finite difference - perturb the part up and down by a
tiny fraction and see how much the dB response moves). Given each
component's assumed tolerance (treated here as a 1-sigma spread),
the standard first-order uncertainty-propagation formula combines
them:

    sigma_dB(f) = sqrt( sum_i ( dDB/dp_i * sigma_i )^2 )   [RSS]

This assumes the component errors are independent and small enough
that the linear (first-order Taylor) approximation holds - true
here since tolerances are a few percent, not tens of percent.
RSS is the "realistic" band: how far a typical real board is
likely to land, since it's very unlikely every component errs in
the same direction simultaneously. We also compute the "worst
case" bound (every component simultaneously at its tolerance limit,
in whichever direction hurts most) as a linear sum instead of RSS -
useful to know as a paranoid upper bound, but expect real boards to
sit well inside it.

Requires: riaa_amp_tf.py in the same directory, numpy, scipy, matplotlib.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from riaa_amp_tf import riaa_amp_tf


def bode_db_deg(G, w):
    """Magnitude (dB) and unwrapped phase (deg) of G at w (rad/s)."""
    _, h = signal.freqs(G.num, G.den, worN=w)
    return 20 * np.log10(np.abs(h)), np.degrees(np.unwrap(np.angle(h)))


# Nominal component values - AS-BUILT standard-value parts, not the
# book's exact numbers (riaa_amp_tf's zero-argument default uses the
# book's exact values; here we use the nearest standard E96/E24 parts
# actually going on the board, per the tolerance/sourcing discussion).
# Reference: Douglas Self's Small Signal Audio Design (2024)
R0 = 1e3;    C0 = 7.96e-6   # electrolytic, loose tolerance regardless of nominal
R1 = 511e3;  C1 = 6.2e-9
R2 = 43.2e3; C2 = 1.8e-9    # update this (and its tolerance below) if you
                            # combine two parts for tighter accuracy here

p0 = np.array([R0, C0, R1, C1, R2, C2])
names = ['R0', 'C0', 'R1', 'C1', 'R2', 'C2']

# Per-component tolerance (fractional, e.g. 0.01 = 1%), treated as 1-sigma
# 1% metal-film resistors, 5% precision film/C0G caps for the RIAA
# timing network, 20% for the electrolytic C0 (typical for electrolytics
# - this breakpoint is the subsonic pole, not one of the core RIAA
# timing constants, so it's far less sensitive to begin with).
tol = np.array([0.01, 0.20, 0.01, 0.05, 0.01, 0.05])  # [R0 C0 R1 C1 R2 C2]

# Frequency sweep
f = np.logspace(np.log10(10), np.log10(100e3), 400)  # 10 Hz to 100 kHz
w = 2 * np.pi * f                                     # rad/s

# Nominal response
G0 = riaa_amp_tf(*p0)
mag0_dB, phase0deg = bode_db_deg(G0, w)

# Sensitivity via central-difference finite differences
n = p0.size
dMagdP = np.zeros((f.size, n))    # d(dB)/d(param), one column per component
dPhasedP = np.zeros((f.size, n))  # d(deg)/d(param)

frac_step = 1e-4  # 0.01% nudge - small enough for a clean numerical
                  # derivative, unrelated to the component's real tolerance

for i in range(n):
    pPlus = p0.copy();  pPlus[i] = p0[i] * (1 + frac_step)
    pMinus = p0.copy(); pMinus[i] = p0[i] * (1 - frac_step)
    h = pPlus[i] - pMinus[i]

    magP, phaseP = bode_db_deg(riaa_amp_tf(*pPlus), w)
    magM, phaseM = bode_db_deg(riaa_amp_tf(*pMinus), w)

    dMagdP[:, i] = (magP - magM) / h
    dPhasedP[:, i] = (phaseP - phaseM) / h

# Propagate uncertainty
sigma = tol * p0  # 1-sigma of each component, in its own native units (ohms/farads)

# RSS (statistical / "realistic spread across real boards") combination.
# Written as a matrix-vector product: for each frequency row, this sums
# (dMagdP_i)^2 * sigma_i^2 over the 6 components in one shot.
sigma_mag_dB = np.sqrt(dMagdP ** 2 @ sigma ** 2)
sigma_phase_deg = np.sqrt(dPhasedP ** 2 @ sigma ** 2)

# Worst-case (linear sum / "every component simultaneously unlucky") bound.
worst_mag_dB = np.abs(dMagdP) @ sigma
worst_phase_deg = np.abs(dPhasedP) @ sigma

# Ideal IEC-amended RIAA target curve, for comparison only
T1 = 3180e-6; T2 = 318e-6; T3 = 75e-6; T4 = 0  # 50Hz, 500Hz, 2122Hz, 20Hz(IEC)
# T4 = 7950e-6
s_jw = 1j * w
Hriaa = (1 + s_jw * T2) / ((1 + s_jw * T1) * (1 + s_jw * T3) * (1 + s_jw * T4))
Hriaa_dB = 20 * np.log10(np.abs(Hriaa))

# Normalize both curves to 0 dB at 1 kHz for a fair shape comparison
idx1k = np.argmin(np.abs(f - 1000))
Hriaa_dB_norm = Hriaa_dB - Hriaa_dB[idx1k]
mag0_dB_norm = mag0_dB - mag0_dB[idx1k]

# --- Plot: magnitude Bode with uncertainty band ---
fig, (ax_mag, ax_ph) = plt.subplots(2, 1, figsize=(9, 7))
fig.canvas.manager.set_window_title(
    'RIAA Amplifier - Bode with Tolerance Error Propagation')

k_sigma = 3  # show +/-3 sigma (~99.7% coverage if errors are ~Gaussian)
band_color = (0.85, 0.90, 1.0)

ax_mag.fill_between(f, mag0_dB_norm - k_sigma * sigma_mag_dB,
                    mag0_dB_norm + k_sigma * sigma_mag_dB,
                    color=band_color, linewidth=0,
                    label=rf'$\pm{k_sigma}\sigma$ (RSS, realistic spread)')
ax_mag.plot(f, mag0_dB_norm + worst_mag_dB, 'r--',
            label='Worst case (all parts unlucky)')
ax_mag.plot(f, mag0_dB_norm - worst_mag_dB, 'r--')
ax_mag.plot(f, Hriaa_dB_norm, 'k:', linewidth=1.5, label='Ideal IEC RIAA target')
ax_mag.plot(f, mag0_dB_norm, 'b-', linewidth=1.5, label='Nominal (as-built values)')

# Literal error bars at standard reference frequencies
f_ref = np.array([20, 50, 100, 500, 1000, 2122, 10000, 20000])
mag_ref = np.interp(f_ref, f, mag0_dB_norm)
sig_ref = np.interp(f_ref, f, sigma_mag_dB)
ax_mag.errorbar(f_ref, mag_ref, yerr=k_sigma * sig_ref, fmt='o',
                color=(0, 0.45, 0), markerfacecolor=(0, 0.45, 0),
                linewidth=1.2,
                label=rf'$\pm{k_sigma}\sigma$ at reference points')

ax_mag.set_xscale('log')
ax_mag.grid(True, which='both')
ax_mag.set_xlabel('Frequency (Hz)')
ax_mag.set_ylabel('Magnitude (dB, normalized to 0 dB @ 1 kHz)')
ax_mag.set_title('RIAA Amplifier Frequency Response with Component Tolerance Uncertainty')
ax_mag.legend(loc='lower left')
ax_mag.set_xlim(10, 100e3)

ax_ph.fill_between(f, phase0deg - k_sigma * sigma_phase_deg,
                   phase0deg + k_sigma * sigma_phase_deg,
                   color=band_color, linewidth=0,
                   label=rf'$\pm{k_sigma}\sigma$ (RSS)')
ax_ph.plot(f, phase0deg, 'b-', linewidth=1.5, label='Nominal')
ax_ph.set_xscale('log')
ax_ph.grid(True, which='both')
ax_ph.set_xlabel('Frequency (Hz)')
ax_ph.set_ylabel('Phase (deg)')
ax_ph.set_title('Phase Response with Tolerance Uncertainty')
ax_ph.legend(loc='best')
ax_ph.set_xlim(10, 100e3)

fig.tight_layout()

# Report which component dominates the error, by frequency band
print('\nDominant error contributor by frequency band '
      '(largest avg |dMag/dp * sigma|):')
bands = [(20, 100), (100, 1000), (1000, 10000), (10000, 20000)]
for lo, hi in bands:
    idxBand = (f >= lo) & (f <= hi)
    contrib = np.mean(np.abs(dMagdP[idxBand, :]) * sigma, axis=0)
    worstIdx = np.argmax(contrib)
    print(f'  {lo:6.0f} - {hi:6.0f} Hz: {names[worstIdx]}  '
          f'(avg contribution {contrib[worstIdx]:.4f} dB)')

plt.show()
