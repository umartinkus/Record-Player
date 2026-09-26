"""CARTRIDGE_PARAM_SWEEP

Python port of cartridge_param_sweep.m.

Maps which combinations of the two UNMEASURED cartridge parameters
are consistent with the AT-VM95E datasheet's published claims.

The Shure/Hunt equivalent circuit needs four parameters. Two we have
data for, two we don't:

  C1  vinyl surface compliance   - property of the record + stylus
                                   contact patch, not the cartridge.
                                   Taken from Pspatial Audio.
  C2  bearing compliance         - AT publishes this (dynamic
                                   compliance @ 100 Hz).
  L1  armature effective mass    - NOT published. Sweep.
  R   bearing viscous damping    - NOT published. Sweep.

Rather than guessing L1 and R, we sweep both and reject any pair
that would contradict something AT actually claims in print. What
survives is the feasible region. This is using published specs to
constrain unpublished parameters.

ACCEPTANCE CRITERIA
  A) HF resonance at or above F_MIN_HF. AT specifies response to
     20 kHz, so the L1/C1 resonance cannot sit below it.
     Depends on L1 only -> a hard ceiling, computed analytically.
  B) Peak tracking force demand <= VTF_MAX at the test velocity.
     Hforce is the downforce the cartridge DEMANDS to stay in the
     groove; AT's nominal VTF is what's available. Demand exceeding
     supply means mistracking.
  C) Velocity response flat within FLAT_TOL_DB over the audio band,
     consistent with a published frequency response spec.

Uses direct complex arithmetic (numpy broadcasting over the whole
L1 x R x f grid) rather than transfer-function objects. Use
cartridge_tf_mech.py for detailed single-point analysis where you want
a TransferFunction object.

NOTE ON UNITS: the direct analogy's own units throughout.
  mH = mg (mass)        uF = u(cm/dyne) (compliance)
  ohm = dyne-s/cm       A = cm/s        V = grams
Conveniently 1 um/mN and 1 u(cm/dyne) are numerically identical
(both 1e-3 m/N), so datasheet compliance drops straight in.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# ------------------------------------------------------------------
#  Known / assumed parameters
# ------------------------------------------------------------------
C1 = 0.0506e-6     # vinyl surface compliance, u(cm/dyne). From Pspatial
                   # Audio, back-calculated from their stated 25 kHz
                   # parallel resonance. Scales roughly as 1/(stylus
                   # contact area) - a Shibata or microline profile in
                   # the same VM95 body would REDUCE this and push the
                   # HF resonance up. Swap here to model a stylus
                   # upgrade.

C2 = 6.5e-6        # bearing compliance, u(cm/dyne) = um/mN.
                   # AT-VM95E published dynamic compliance @ 100 Hz.
                   # NOTE: this is the 100 Hz figure. Viscoelastic
                   # compliance rises at lower frequency - the datasheet
                   # static figure is larger. See the Zener extension in
                   # the notes at the bottom of this file.

# ------------------------------------------------------------------
#  Acceptance criteria - EDIT THESE against the datasheet
# ------------------------------------------------------------------
F_MIN_HF = 20e3      # Hz. Lowest acceptable L1/C1 resonance.
VTF_MAX = 2.0        # grams. AT-VM95E nominal vertical tracking force.
V_TEST = 5.0         # cm/s. Test velocity for the force check. 5 cm/s
                     # is the reference level; real test records reach
                     # 15-30 cm/s at HF, so treat this as a baseline
                     # and re-run at higher values for a stress test.
FLAT_TOL_DB = 2.0    # dB. Allowed p-p ripple in Hvel, 20 Hz - 20 kHz.

# ------------------------------------------------------------------
#  Sweep ranges
# ------------------------------------------------------------------
L1_vec = np.linspace(0.2e-3, 1.4e-3, 120)  # mg. Typical MM tip mass is
                                           # 0.3-1.0 mg for an aluminium
                                           # cantilever.
R_vec = np.linspace(2, 120, 120)           # dyne-s/cm.

f = np.logspace(np.log10(20), np.log10(20e3), 300)
w = 2 * np.pi * f
jw = 1j * w

# ------------------------------------------------------------------
#  Analytical bound from criterion A (no sweep needed)
# ------------------------------------------------------------------
L1_ceiling = 1 / (C1 * (2 * np.pi * F_MIN_HF) ** 2)
print(f'Criterion A gives a hard ceiling: L1 <= {L1_ceiling * 1e3:.3f} mg')
print(f'  (any heavier armature puts the HF resonance below '
      f'{F_MIN_HF / 1e3:.0f} kHz)\n')

# Critical damping of the numerator - below this R the force notch
# becomes a real complex-conjugate notch rather than two real zeros.
print('For reference, R_crit = 2*sqrt(L1/C2):')  # series RLC resonance
print(f'  at L1 = 0.8 mg -> R_crit = {2 * np.sqrt(0.8e-3 / C2):.1f} dyne-s/cm\n')

# ------------------------------------------------------------------
#  Sweep
# ------------------------------------------------------------------
# Grid axes: [L1 index, R index, frequency index]
L1g = L1_vec[:, None, None]
Rg = R_vec[None, :, None]

Zc1 = 1 / (jw * C1)                       # independent of the swept params
Zbranch = jw * L1g + 1 / (jw * C2) + Rg

Hvel = Zc1 / (Zc1 + Zbranch)              # tf from Ii to Io
Hforce = (Zc1 * Zbranch) / (Zc1 + Zbranch)  # Ii to Vab

# HF resonance, Hz - every column is identical (depends on L1 only)
f_hf = np.broadcast_to(1 / (2 * np.pi * np.sqrt(L1_vec * C1))[:, None],
                       (L1_vec.size, R_vec.size))

# Hforce is grams per (cm/s); scale by the test velocity.
peak_force = np.max(np.abs(Hforce), axis=2) * V_TEST  # grams

magdB = 20 * np.log10(np.abs(Hvel))
ripple_dB = magdB.max(axis=2) - magdB.min(axis=2)     # p-p flatness, dB

# ------------------------------------------------------------------
#  Feasibility
# ------------------------------------------------------------------
passA = f_hf >= F_MIN_HF
passB = peak_force <= VTF_MAX
passC = ripple_dB <= FLAT_TOL_DB
feasible = passA & passB & passC

print(f'Feasible fraction of swept space: {100 * feasible.mean():.1f}%')
if not feasible.any():
    print('  NOTHING PASSES. Either a criterion is too tight, or the\n'
          '  sweep range misses the feasible region, or the model is\n'
          '  wrong. Check the individual pass maps below before\n'
          '  loosening anything.')

# ------------------------------------------------------------------
#  Plots
# ------------------------------------------------------------------
RR, LL = np.meshgrid(R_vec, L1_vec * 1e3)  # L1 in mg for readability

fig, axs = plt.subplots(2, 2, figsize=(10, 7.5))
fig.canvas.manager.set_window_title('Cartridge parameter feasibility')


def _contour_panel(ax, Z, level, title, Z_limit=None):
    cf = ax.contourf(RR, LL, Z, 20)
    ax.contour(RR, LL, Z if Z_limit is None else Z_limit, [level],
               colors='w', linewidths=2)
    fig.colorbar(cf, ax=ax)
    ax.set_xlabel('R (dyne-s/cm)')
    ax.set_ylabel('L1 (mg)')
    ax.set_title(title)


# f_hf varies along L1 only, so its contour levels are horizontal lines.
_contour_panel(axs[0, 0], f_hf / 1e3, F_MIN_HF,
               'A: HF resonance (kHz), white = limit', Z_limit=f_hf)
_contour_panel(axs[0, 1], peak_force, VTF_MAX,
               f'B: peak force demand (g) at {V_TEST:.0f} cm/s')
_contour_panel(axs[1, 0], ripple_dB, FLAT_TOL_DB,
               'C: response ripple (dB p-p)')

ax = axs[1, 1]
ax.pcolormesh(R_vec, L1_vec * 1e3, feasible.astype(float), shading='auto',
              cmap=ListedColormap([[0.9, 0.9, 0.9], [0.2, 0.6, 0.4]]),
              vmin=0, vmax=1)
ax.set_xlabel('R (dyne-s/cm)')
ax.set_ylabel('L1 (mg)')
ax.set_title('Feasible region (all three criteria)')

fig.tight_layout()

# ------------------------------------------------------------------
#  Report the centroid of the feasible region as a working estimate
# ------------------------------------------------------------------
if feasible.any():
    iF, jF = np.nonzero(feasible)
    L1_est = L1_vec[iF].mean()
    R_est = R_vec[jF].mean()
    print('\nCentroid of feasible region (use as working values):')
    print(f'  L1 = {L1_est * 1e3:.3f} mg')
    print(f'  R  = {R_est:.1f} dyne-s/cm')
    print('\nRanges spanned by the feasible region:')
    print(f'  L1: {L1_vec[iF].min() * 1e3:.3f} to {L1_vec[iF].max() * 1e3:.3f} mg')
    print(f'  R : {R_vec[jF].min():.1f} to {R_vec[jF].max():.1f} dyne-s/cm')
    print('\nThe centroid is a convenience, not a physical estimate -\n'
          'it has no meaning if the feasible region is L-shaped or\n'
          'disconnected. Look at the map before trusting it.')

plt.show()

# ------------------------------------------------------------------
#  NOTES
# ------------------------------------------------------------------
# Once you have the datasheet static compliance, replace the fixed C2
# with a standard-linear-solid (Zener) bearing, which stays rational:
#
#     Ceff(s) = Cinf + (C0 - Cinf)/(1 + s*tau)
#     Zbear(s) = 1/(s*Ceff(s))
#
# and drop R, since the relaxation IS the damping - keeping both
# double-counts the same physical mechanism. C0 comes from the static
# spec, Cinf from the 100 Hz dynamic spec, leaving tau as the single
# free parameter to sweep. That turns this 2D sweep into a 2D sweep
# over (L1, tau) with better physical grounding.
#
# For Monte Carlo instead of a sweep: sample L1 and R from distributions
# over the feasible region found here, propagate through the full
# cascade (mech * k * elec), and look at the spread in end-to-end
# response. Do that AFTER this sweep, so you know which parameters
# actually move the output - you may find only one of them matters and
# the Monte Carlo collapses to 1D.
