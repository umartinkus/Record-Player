"""Shure/Hunt phono cartridge equivalent circuit.

Python port of cartridge_tf_mech.m.

Shure/Hunt phono cartridge equivalent circuit (Pspatial Audio, after
Anderson et al., JAES 1966). Direct analogy: force~V, velocity~I,
mass~L, compliance~C, damping~R.

Units are the analogy's own: mH=mg, uF=u(cm/dyne), ohm=dyne-s/cm,
A=cm/s, V=grams. Work numerically in electrical units and results
come out in the mechanical ones.

  Hvel   = Io/Ii, stylus velocity / recorded velocity -> freq response
  Hforce = Vab/Ii, tracking force / recorded velocity [grams/(cm/s)]

Topology: current source Ii at node A; C1 (vinyl compliance) shunts A
to ground; series branch A -> L1 -> C2 -> R -> ground. R sits in series
with C2 because the elastomer bearing is viscoelastic and the two are
physically inseparable.
"""

from scipy import signal


def cartridge_tf_mech(L1_mH=0.80e-3, C1_uF=0.0506e-6, C2_uF=25e-6, R_ohm=40):
    """Return (Hvel, Hforce) as scipy.signal.TransferFunction objects.

    Call with no arguments for the book defaults, or pass all four.

    MATLAB builds these with tf('s') algebra and minreal(); here the
    reduced forms are written out directly:

        Zbranch = s*L1 + 1/(s*C2) + R = (L1*C2*s^2 + R*C2*s + 1) / (s*C2)
        Zc1     = 1/(s*C1)

        Hvel   = Zc1 / (Zc1 + Zbranch) = 1 / (1 + s*C1*Zbranch)
               = C2 / (C1*L1*C2*s^2 + C1*R*C2*s + (C1 + C2))

        Hforce = Zc1*Zbranch / (Zc1 + Zbranch) = Zbranch * Hvel
               = (L1*C2*s^2 + R*C2*s + 1)
                 / (s * (C1*L1*C2*s^2 + C1*R*C2*s + (C1 + C2)))
    """
    L1, C1, C2, R = L1_mH, C1_uF, C2_uF, R_ohm

    den = [C1 * L1 * C2, C1 * R * C2, C1 + C2]
    Hvel = signal.TransferFunction([C2], den)
    Hforce = signal.TransferFunction([L1 * C2, R * C2, 1.0], den + [0.0])
    return Hvel, Hforce
