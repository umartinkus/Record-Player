"""Laplace-domain transfer function of the single-stage, non-inverting
active RIAA amplifier (Douglas Self, "Small Signal Audio Design").

Python port of riaa_amp_tf.m.

This models the OP-AMP GAIN STAGE ONLY: input at the op-amp's
non-inverting (+) terminal, output at the op-amp output. It does NOT
include the cartridge's source impedance, 47k/pF loading network, or
cable capacitance - those belong in a separate front-end block when you
assemble the full system model.

  G = riaa_amp_tf() returns the transfer function using the exact
      reference component values from Self's book:
          R0 = 1e3        C0 = 7.96e-6
          R1 = 513.5e3    C1 = 6.193e-9
          R2 = 42.86e3    C2 = 1.75e-9

  G = riaa_amp_tf(R0, C0, R1, C1, R2, C2) uses custom component values
      (SI units throughout: ohms, farads). Call this with perturbed
      values for tolerance sweeps, Monte Carlo analysis, or to model
      the as-built board using real (standard-value) parts instead
      of the book's exact numbers.

  Circuit (non-inverting op-amp topology):

      Gain(s) = 1 + Zf(s) / Zg(s)

      Zg(s) = R0 + 1/(s*C0)                      ground leg, from the
                                                   inverting (-) input
      Zf(s) = R1/(1+s*R1*C1) + R2/(1+s*R2*C2)     feedback network,
                                                   output to (-) input

  With a finite open-loop gain A0 the closed-loop gain becomes

      G(s) = A0 / (1 + A0 * Zg/(Zg + Zf))

  R0*C0 sets the subsonic/DC-gain-unity time constant (~20 Hz, the
  IEC amendment pole - it's what keeps DC gain at unity and rolls off
  turntable rumble, built into the network rather than a separate
  stage). R1*C1 = 3180us (50 Hz pole). R2*C2 = 75us (2122 Hz pole).
  The 500 Hz breakpoint is not a direct RC product - it emerges from
  the interaction of the two feedback branches (see Lipshitz, "On
  RIAA Equalization Networks", JAES 1979, for the full derivation).

  Returns G as a scipy.signal.TransferFunction (Laplace/s domain).

  Example:
      from scipy import signal
      import numpy as np

      G = riaa_amp_tf()
      w, mag_dB, phase_deg = signal.bode(G)

      # gain at 1 kHz, in dB
      _, h = signal.freqs(G.num, G.den, worN=[2*np.pi*1000])
      g1k_dB = 20*np.log10(abs(h[0]))

      # swap in as-built standard-value parts instead of the book's
      # exact values
      G_asbuilt = riaa_amp_tf(1e3, 8.2e-6, 511e3, 6.2e-9, 43.2e3, 1.8e-9)

  See also: riaa_bode_tolerance_analysis.py (separate analysis script
  that calls this function to build a Bode plot with a component-
  tolerance uncertainty band via error propagation).
"""

import numpy as np
from scipy import signal


def riaa_amp_tf(R0=1e3, C0=7.96e-6, R1=513.5e3, C1=6.193e-9,
                R2=42.86e3, C2=1.75e-9):
    A0 = 2e6

    # Clear every impedance's denominator by multiplying through by
    # s*C0*(1 + s*R1*C1)*(1 + s*R2*C2). With
    #     Ng   = (s*R0*C0 + 1)(1 + s*R1*C1)(1 + s*R2*C2)      (Zg term)
    #     Nsum = Ng + s*C0*R1*(1 + s*R2*C2) + s*C0*R2*(1 + s*R1*C1)
    # the closed-loop gain is
    #     G = A0*Nsum / (Nsum + A0*Ng)
    # which is the already-minimal form minreal() produces in MATLAB.
    t1 = [R1 * C1, 1.0]
    t2 = [R2 * C2, 1.0]
    Ng = np.polymul(np.polymul([R0 * C0, 1.0], t1), t2)
    Nsum = np.polyadd(Ng, np.polyadd(np.polymul([C0 * R1, 0.0], t2),
                                     np.polymul([C0 * R2, 0.0], t1)))

    num = A0 * Nsum
    den = np.polyadd(Nsum, A0 * Ng)
    return signal.TransferFunction(num, den)
