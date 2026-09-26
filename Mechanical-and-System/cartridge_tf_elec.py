"""VM95E electrical network: EMF -> preamp input.

Python port of cartridge_tf_elec.m.

Generator is flat (Faraday: EMF proportional to velocity); all the
frequency shaping lives in Lcoil/Rcoil driving Rload || Ctot.

  Hel - EMF at generator -> voltage at preamp input (dimensionless)
  k   - sensitivity, V per (m/s) of stylus velocity

Lcoil - coil inductance, Rcoil - coil resistance
Ctot - capacitance of tone arm interconnect
Defaults: Lcoil=0.55, Rcoil=485, Rload=47e3, Ctot=100e-12..200e-12
Ctot is EVERYTHING: tonearm wiring + interconnect + preamp input.
"""

from scipy import signal


def cartridge_tf_elec(Lcoil_H=550e-3, Rcoil_ohm=485, Rload_ohm=47000,
                      Ctot_F=100e-12):
    """Return (Hel, k). Hel is a scipy.signal.TransferFunction.

    Call with no arguments for the book defaults, or pass all four.

        Zp  = Rload / (1 + s*Rload*Ctot)            load in parallel with C
        Hel = Zp / (Zp + s*Lcoil + Rcoil)
            = Rload / (Lcoil*Rload*Ctot*s^2
                       + (Lcoil + Rcoil*Rload*Ctot)*s + (Rcoil + Rload))
    """
    L, Rc, Rl, C = Lcoil_H, Rcoil_ohm, Rload_ohm, Ctot_F

    Hel = signal.TransferFunction([Rl], [L * Rl * C, L + Rc * Rl * C, Rc + Rl])
    k = 0.08  # 4mV / 5cm/s -> V/(m/s)
    return Hel, k
