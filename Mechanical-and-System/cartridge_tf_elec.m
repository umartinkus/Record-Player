function [Hel, k] = cartridge_tf_elec(Lcoil_H, Rcoil_ohm, Rload_ohm, Ctot_F)
% VM95E electrical network: EMF -> preamp input.
% Generator is flat (Faraday: EMF proportional to velocity); all the
% frequency shaping lives in Lcoil/Rcoil driving Rload || Ctot.
%
%   Hel - EMF at generator -> voltage at preamp input (dimensionless)
%   k   - sensitivity, V per (m/s) of stylus velocity
%
% Lcoil - coil inductance, Rcoil - coil resistance
% Ctot - capacitance of tone arm interconnect
% Defaults: Lcoil=0.55, Rcoil=485, Rload=47e3, Ctot=100e-12..200e-12
% Ctot is EVERYTHING: tonearm wiring + interconnect + preamp input.
    if nargin == 0
        Rload_ohm = 47000; Lcoil_H = 550e-3;
        Ctot_F = 100e-12; Rcoil_ohm = 485;

    elseif nargin ~= 4
        error('cartridge_tf_mech:nargin', ...
            'Call with zero arguments (book defaults) or all args.');
    end
    s  = tf('s');
    Zp = Rload_ohm / (1 + s*Rload_ohm*Ctot_F);          % load in parallel with C
    Hel = minreal(Zp / (Zp + s*Lcoil_H + Rcoil_ohm));
    k   = 0.08;                                % 4mV / 5cm/s -> V/(m/s)
end