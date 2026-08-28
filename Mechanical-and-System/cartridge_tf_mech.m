function [Hvel, Hforce] = cartridge_tf_mech(L1_mH, C1_uF, C2_uF, R_ohm)
% Shure/Hunt phono cartridge equivalent circuit (Pspatial Audio, after
% Anderson et al., JAES 1966). Direct analogy: force~V, velocity~I,
% mass~L, compliance~C, damping~R.
%
% Units are the analogy's own: mH=mg, uF=u(cm/dyne), ohm=dyne-s/cm,
% A=cm/s, V=grams. Work numerically in electrical units and results
% come out in the mechanical ones.
%
%   Hvel   = Io/Ii, stylus velocity / recorded velocity -> freq response
%   Hforce = Vab/Ii, tracking force / recorded velocity [grams/(cm/s)]
%
% Topology: current source Ii at node A; C1 (vinyl compliance) shunts A
% to ground; series branch A -> L1 -> C2 -> R -> ground. R sits in series
% with C2 because the elastomer bearing is viscoelastic and the two are
% physically inseparable.
    if nargin == 0
        L1_mH = 0.80e-3; C2_uF = 25e-6;
        C1_uF = 0.0506e-6; R_ohm = 40;
    elseif nargin ~= 4
        error('cartridge_tf_mech:nargin', ...
            'Call with zero arguments (book defaults) or all args.');
    end
    s = tf('s');
    Zbranch = s*L1_mH + 1/(s*C2_uF) + R_ohm; % armature + bearing + damping
    Zc1     = 1/(s*C1_uF); % vinyl surface compliance

    Hvel   = minreal(Zc1 / (Zc1 + Zbranch));
    Hforce = minreal((Zc1*Zbranch) / (Zc1 + Zbranch));
end