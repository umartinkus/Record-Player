%CARTRIDGE_PARAM_SWEEP
%   Maps which combinations of the two UNMEASURED cartridge parameters
%   are consistent with the AT-VM95E datasheet's published claims.
%
%   The Shure/Hunt equivalent circuit needs four parameters. Two we have
%   data for, two we don't:
%
%     C1  vinyl surface compliance   - property of the record + stylus
%                                      contact patch, not the cartridge.
%                                      Taken from Pspatial Audio.
%     C2  bearing compliance         - AT publishes this (dynamic
%                                      compliance @ 100 Hz).
%     L1  armature effective mass    - NOT published. Sweep.
%     R   bearing viscous damping    - NOT published. Sweep.
%
%   Rather than guessing L1 and R, we sweep both and reject any pair
%   that would contradict something AT actually claims in print. What
%   survives is the feasible region. This is using published specs to
%   constrain unpublished parameters.
%
%   ACCEPTANCE CRITERIA
%     A) HF resonance at or above F_MIN_HF. AT specifies response to
%        20 kHz, so the L1/C1 resonance cannot sit below it.
%        Depends on L1 only -> a hard ceiling, computed analytically.
%     B) Peak tracking force demand <= VTF_MAX at the test velocity.
%        Hforce is the downforce the cartridge DEMANDS to stay in the
%        groove; AT's nominal VTF is what's available. Demand exceeding
%        supply means mistracking.
%     C) Velocity response flat within FLAT_TOL_DB over the audio band,
%        consistent with a published frequency response spec.
%
%   Uses direct complex arithmetic rather than tf objects - roughly two
%   orders of magnitude faster inside a nested loop, and it drops the
%   Control System Toolbox dependency for this script. Use cartridge_tf.m
%   for detailed single-point analysis where you want a tf object.
%
%   NOTE ON UNITS: the direct analogy's own units throughout.
%     mH = mg (mass)        uF = u(cm/dyne) (compliance)
%     ohm = dyne-s/cm       A = cm/s        V = grams
%   Conveniently 1 um/mN and 1 u(cm/dyne) are numerically identical
%   (both 1e-3 m/N), so datasheet compliance drops straight in.

clear; clc; close all;

%% ------------------------------------------------------------------
%  Known / assumed parameters
%  ------------------------------------------------------------------
C1 = 0.0506e-6;    % vinyl surface compliance, u(cm/dyne). From Pspatial
                   % Audio, back-calculated from their stated 25 kHz
                   % parallel resonance. Scales roughly as 1/(stylus
                   % contact area) - a Shibata or microline profile in
                   % the same VM95 body would REDUCE this and push the
                   % HF resonance up. Swap here to model a stylus
                   % upgrade.

C2 = 6.5e-6;       % bearing compliance, u(cm/dyne) = um/mN.
                   % AT-VM95E published dynamic compliance @ 100 Hz.
                   % NOTE: this is the 100 Hz figure. Viscoelastic
                   % compliance rises at lower frequency - the datasheet
                   % static figure is larger. See the Zener extension in
                   % the notes at the bottom of this file.

%% ------------------------------------------------------------------
%  Acceptance criteria - EDIT THESE against the datasheet
%  ------------------------------------------------------------------
F_MIN_HF    = 20e3;   % Hz. Lowest acceptable L1/C1 resonance.
VTF_MAX     = 2.0;    % grams. AT-VM95E nominal vertical tracking force.
V_TEST      = 5.0;    % cm/s. Test velocity for the force check. 5 cm/s
                      % is the reference level; real test records reach
                      % 15-30 cm/s at HF, so treat this as a baseline
                      % and re-run at higher values for a stress test.
FLAT_TOL_DB = 2.0;    % dB. Allowed p-p ripple in Hvel, 20 Hz - 20 kHz.

%% ------------------------------------------------------------------
%  Sweep ranges
%  ------------------------------------------------------------------
L1_vec = linspace(0.2e-3, 1.4e-3, 120);   % mg. Typical MM tip mass is
                                          % 0.3-1.0 mg for an aluminium
                                          % cantilever.
R_vec  = linspace(2, 120, 120);           % dyne-s/cm.

f = logspace(log10(20), log10(20e3), 300);
w = 2*pi*f;
jw = 1i*w;

%% ------------------------------------------------------------------
%  Analytical bound from criterion A (no sweep needed)
%  ------------------------------------------------------------------
L1_ceiling = 1 / (C1 * (2*pi*F_MIN_HF)^2);
fprintf('Criterion A gives a hard ceiling: L1 <= %.3f mg\n', L1_ceiling*1e3);
fprintf('  (any heavier armature puts the HF resonance below %.0f kHz)\n\n', ...
        F_MIN_HF/1e3);

% Critical damping of the numerator - below this R the force notch
% becomes a real complex-conjugate notch rather than two real zeros.
fprintf('For reference, R_crit = 2*sqrt(L1/C2):\n'); % series RLC resonance
fprintf('  at L1 = 0.8 mg -> R_crit = %.1f dyne-s/cm\n\n', ...
        2*sqrt(0.8e-3/C2));

%% ------------------------------------------------------------------
%  Sweep
%  ------------------------------------------------------------------
nL = numel(L1_vec);  nR = numel(R_vec);

f_hf      = zeros(nL, nR);   % HF resonance, Hz
peak_force = zeros(nL, nR);  % max demanded tracking force, grams
ripple_dB  = zeros(nL, nR);  % p-p flatness of Hvel, dB

Zc1 = 1 ./ (jw * C1);        % independent of the swept params

for i = 1:nL
    L1 = L1_vec(i);
    for j = 1:nR
        R = R_vec(j);

        Zbranch = jw*L1 + 1./(jw*C2) + R;

        Hvel   = Zc1 ./ (Zc1 + Zbranch); % tf from Ii to Io
        Hforce = (Zc1 .* Zbranch) ./ (Zc1 + Zbranch); % Ii to Vab

        f_hf(i,j) = 1 / (2*pi*sqrt(L1*C1)); % every col is linearly dep

        % Hforce is grams per (cm/s); scale by the test velocity.
        peak_force(i,j) = max(abs(Hforce)) * V_TEST;

        magdB = 20*log10(abs(Hvel));
        ripple_dB(i,j) = max(magdB) - min(magdB);
    end
end

%% ------------------------------------------------------------------
%  Feasibility
%  ------------------------------------------------------------------
passA = f_hf       >= F_MIN_HF;
passB = peak_force <= VTF_MAX;
passC = ripple_dB  <= FLAT_TOL_DB;
feasible = passA & passB & passC;

fprintf('Feasible fraction of swept space: %.1f%%\n', 100*mean(feasible(:)));
if ~any(feasible(:))
    fprintf(['  NOTHING PASSES. Either a criterion is too tight, or the\n' ...
             '  sweep range misses the feasible region, or the model is\n' ...
             '  wrong. Check the individual pass maps below before\n' ...
             '  loosening anything.\n']);
end

%% ------------------------------------------------------------------
%  Plots
%  ------------------------------------------------------------------
[RR, LL] = meshgrid(R_vec, L1_vec*1e3);   % L1 in mg for readability

figure('Name','Cartridge parameter feasibility','Position',[100 100 1000 750]);

subplot(2,2,1);
contourf(RR, LL, f_hf/1e3, 20, 'LineColor','none'); hold on;
contour(RR, LL, f_hf, [F_MIN_HF F_MIN_HF], 'w-', 'LineWidth', 2);
colorbar; xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title('A: HF resonance (kHz), white = limit');

subplot(2,2,2);
contourf(RR, LL, peak_force, 20, 'LineColor','none'); hold on;
contour(RR, LL, peak_force, [VTF_MAX VTF_MAX], 'w-', 'LineWidth', 2);
colorbar; xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title(sprintf('B: peak force demand (g) at %.0f cm/s', V_TEST));

subplot(2,2,3);
contourf(RR, LL, ripple_dB, 20, 'LineColor','none'); hold on;
contour(RR, LL, ripple_dB, [FLAT_TOL_DB FLAT_TOL_DB], 'w-', 'LineWidth', 2);
colorbar; xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title('C: response ripple (dB p-p)');

subplot(2,2,4);
imagesc(R_vec, L1_vec*1e3, double(feasible)); set(gca,'YDir','normal');
colormap(gca, [0.9 0.9 0.9; 0.2 0.6 0.4]);
xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title('Feasible region (all three criteria)');

%% ------------------------------------------------------------------
%  Report the centroid of the feasible region as a working estimate
%  ------------------------------------------------------------------
if any(feasible(:))
    [iF, jF] = find(feasible);
    L1_est = mean(L1_vec(iF));
    R_est  = mean(R_vec(jF));
    fprintf('\nCentroid of feasible region (use as working values):\n');
    fprintf('  L1 = %.3f mg\n', L1_est*1e3);
    fprintf('  R  = %.1f dyne-s/cm\n', R_est);
    fprintf('\nRanges spanned by the feasible region:\n');
    fprintf('  L1: %.3f to %.3f mg\n', min(L1_vec(iF))*1e3, max(L1_vec(iF))*1e3);
    fprintf('  R : %.1f to %.1f dyne-s/cm\n', min(R_vec(jF)), max(R_vec(jF)));
    fprintf(['\nThe centroid is a convenience, not a physical estimate -\n' ...
             'it has no meaning if the feasible region is L-shaped or\n' ...
             'disconnected. Look at the map before trusting it.\n']);
end

%% ------------------------------------------------------------------
%  NOTES
%  ------------------------------------------------------------------
% Once you have the datasheet static compliance, replace the fixed C2
% with a standard-linear-solid (Zener) bearing, which stays rational:
%
%     Ceff(s) = Cinf + (C0 - Cinf)/(1 + s*tau)
%     Zbear(s) = 1/(s*Ceff(s))
%
% and drop R, since the relaxation IS the damping - keeping both
% double-counts the same physical mechanism. C0 comes from the static
% spec, Cinf from the 100 Hz dynamic spec, leaving tau as the single
% free parameter to sweep. That turns this 2D sweep into a 2D sweep
% over (L1, tau) with better physical grounding.
%
% For Monte Carlo instead of a sweep: sample L1 and R from distributions
% over the feasible region found here, propagate through the full
% cascade (mech * k * elec), and look at the spread in end-to-end
% response. Do that AFTER this sweep, so you know which parameters
% actually move the output - you may find only one of them matters and
% the Monte Carlo collapses to 1D.
