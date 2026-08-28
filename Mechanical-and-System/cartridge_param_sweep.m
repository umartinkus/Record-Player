%CARTRIDGE_PARAM_SWEEP_V2
%   Maps which combinations of the two UNMEASURED cartridge parameters
%   (L1, R) are consistent with the AT-VM95E datasheet's published claims.
%
%   ==================================================================
%   CHANGES FROM V1 - three bugs, all of which made the sweep useless
%   ==================================================================
%
%   1. UNITS. Hforce comes out in DYNE per (cm/s), not gram-force.
%      Compliance in cm/dyne and velocity in cm/s make the analogy's
%      impedance dyne-s/cm. V1 compared dynes against a 2.0 gram-force
%      limit - a factor of 980.7 too strict. Fixed by DYNE_PER_GF.
%
%   2. TEST SIGNAL. V1 used a constant 5 cm/s across the whole band.
%      At 20 Hz that implies 5/(2*pi*20) = 398 um of displacement.
%      Real groove modulation peaks near 50 um, so V1 was demanding
%      about 8x the maximum a cutting lathe can physically cut. Force
%      scales as displacement/compliance, so this demanded ~6.2 gf at
%      20 Hz REGARDLESS of L1 and R - criterion B failed everywhere for
%      a reason unrelated to the swept parameters. Fixed by the
%      three-limit velocity envelope below.
%
%   3. PLOTTING. Peak force spans orders of magnitude (Q ~ 63 at R = 2),
%      so contourf put all its levels in one corner. Fixed with log
%      scaling on the force map, explicit clim() on the feasibility map,
%      and margin plots that show HOW MUCH each criterion passes or
%      fails by rather than just whether it does.
%
%   ==================================================================
%   THE VELOCITY ENVELOPE
%   ==================================================================
%   Recorded velocity is bounded by three independent physical limits:
%
%       v(f) = min( 2*pi*f*A_MAX ,  V_MAX ,  a_max/(2*pi*f) )
%              \_____________/      \___/    \____________/
%               displacement-       velocity- acceleration-
%               limited (groove      limited  limited (cutter
%               pitch), LF           midband  and stylus), HF
%
%   With the defaults below the corners land at ~477 Hz and ~10.4 kHz.
%
%   Consequence worth understanding: in the constant-amplitude region
%   the demanded force is A_MAX/C2 - CONSTANT, and independent of both
%   swept parameters. So the LF end stops discriminating and criterion B
%   starts biting where it should, at the HF resonance, where it
%   discriminates on R through the resonance Q.
%
%   ==================================================================
%   UNITS (the direct analogy's own)
%     mH = mg (mass)          uF = u(cm/dyne) (compliance)
%     ohm = dyne-s/cm         A = cm/s        V = dyne
%   1 um/mN and 1 u(cm/dyne) are numerically identical (both 1e-3 m/N),
%   so datasheet compliance drops straight in with no conversion.
%   ==================================================================

clear; clc; close all;

DYNE_PER_GF = 980.665;   % gravitational acceleration, cm/s^2

%% ------------------------------------------------------------------
%  Known parameters (NOT swept)
%  ------------------------------------------------------------------
C1 = 0.0506e-6;    % vinyl surface compliance, u(cm/dyne).
                   % Property of the record + stylus contact patch, not
                   % the cartridge. Back-calculated from Pspatial's
                   % 25 kHz parallel resonance. Scales roughly as
                   % 1/(contact area) - reduce it to model a Shibata or
                   % microline stylus in the same VM95 body.

C2 = 6.5e-6;       % bearing compliance, u(cm/dyne) = um/mN.
                   % AT-VM95E published DYNAMIC compliance @ 100 Hz.
                   % This is a known datasheet value and is never swept.
                   % (See the Zener note at the bottom for how the
                   % datasheet STATIC compliance enters.)

%% ------------------------------------------------------------------
%  Velocity envelope - what a record can physically carry
%  ------------------------------------------------------------------
A_MAX = 50e-4;     % cm. Peak groove displacement, LF limit. ~50 um is
                   % a loud LP; variable-pitch cutting can exceed it
                   % briefly. This one parameter sets the entire LF
                   % force demand (= A_MAX/C2), so it matters.
V_MAX = 15;        % cm/s. Peak midband velocity. 5 cm/s is the
                   % reference level; 15 is a loud passage.
A_ACC = 1000 * DYNE_PER_GF;   % cm/s^2. Peak acceleration, HF limit.
                   % ~1000 g is the usual figure quoted for tracking-
                   % test records.

%% ------------------------------------------------------------------
%  Acceptance criteria
%  ------------------------------------------------------------------
F_MIN_HF    = 20e3;   % Hz. Lowest acceptable L1/C1 resonance, since AT
                      % specifies response to 20 kHz.
VTF_MAX     = 2.0;    % gram-force. AT-VM95E nominal tracking force.
FLAT_TOL_DB = 5.0;    % dB p-p, 20 Hz - 20 kHz.
                      % *** WEAKEST CRITERION - THIS IS AN ASSUMPTION. ***
                      % AT publishes a response RANGE (20-22 kHz) but no
                      % tolerance, so there is no published number to
                      % test against. Real MM cartridges show +2 to +5 dB
                      % of rise near the tip-mass resonance and that is
                      % normal, not a defect. V1 used 2 dB, which rejects
                      % essentially every physically real cartridge.
                      % Treat results from this criterion with suspicion.

%% ------------------------------------------------------------------
%  Sweep ranges
%  ------------------------------------------------------------------
L1_vec = linspace(0.2e-3, 1.4e-3, 140);   % mg. Typical MM tip mass is
                                          % 0.3-1.0 mg (aluminium
                                          % cantilever).
R_vec  = linspace(5, 200, 140);           % dyne-s/cm. Extended upward
                                          % from V1 - the feasible
                                          % region turns out to sit
                                          % higher than the article's
                                          % R = 40.

f  = logspace(log10(20), log10(20e3), 400);
w  = 2*pi*f;
jw = 1i*w;

%% ------------------------------------------------------------------
%  Build the velocity envelope
%  ------------------------------------------------------------------
v_env = min( min(w*A_MAX, V_MAX), A_ACC./w );   % cm/s vs frequency

f_corner_lo = V_MAX / (2*pi*A_MAX);
f_corner_hi = A_ACC / (2*pi*V_MAX);
fprintf('Velocity envelope corners: %.0f Hz and %.1f kHz\n', ...
        f_corner_lo, f_corner_hi/1e3);
fprintf('LF force demand (displacement-limited region):\n');
fprintf('  A_MAX/C2 = %.0f dyne = %.2f gf  (independent of L1 and R)\n', ...
        A_MAX/C2, A_MAX/C2/DYNE_PER_GF);
if A_MAX/C2/DYNE_PER_GF > VTF_MAX
    fprintf(['  *** This ALONE exceeds VTF_MAX. Nothing can pass.\n' ...
             '      Either A_MAX is too large or the cartridge genuinely\n' ...
             '      cannot track this record at this VTF. ***\n']);
end
fprintf('\n');

%% ------------------------------------------------------------------
%  Analytical bounds (no sweep needed)
%  ------------------------------------------------------------------
L1_ceiling = 1 / (C1 * (2*pi*F_MIN_HF)^2);
fprintf('Criterion A hard ceiling: L1 <= %.3f mg\n', L1_ceiling*1e3);
fprintf('  (heavier armature puts the HF resonance below %.0f kHz)\n\n', ...
        F_MIN_HF/1e3);

fprintf('R_crit = 2*sqrt(L1/C2), below which the force zeros go complex\n');
fprintf('and a real notch appears:\n');
for L1t = [0.4e-3 0.8e-3 1.2e-3]
    fprintf('  L1 = %.1f mg -> R_crit = %.1f\n', L1t*1e3, 2*sqrt(L1t/C2));
end
fprintf('\n');

%% ------------------------------------------------------------------
%  Sweep
%  ------------------------------------------------------------------
nL = numel(L1_vec);  nR = numel(R_vec);

f_hf       = zeros(nL, nR);   % HF resonance, Hz
peak_force = zeros(nL, nR);   % max demanded force, GRAM-FORCE
ripple_dB  = zeros(nL, nR);   % p-p flatness of Hvel, dB
f_at_peak  = zeros(nL, nR);   % where the force peak occurs, Hz

Zc1 = 1 ./ (jw * C1);         % independent of swept params

for i = 1:nL
    L1 = L1_vec(i);
    for j = 1:nR
        R = R_vec(j);

        Zbranch = jw*L1 + 1./(jw*C2) + R;

        Hvel   = Zc1 ./ (Zc1 + Zbranch);
        Hforce = (Zc1 .* Zbranch) ./ (Zc1 + Zbranch);

        f_hf(i,j) = 1 / (2*pi*sqrt(L1*C1));

        % Force demanded by the envelope, converted to gram-force.
        force_gf = abs(Hforce) .* v_env / DYNE_PER_GF;
        [peak_force(i,j), idx] = max(force_gf);
        f_at_peak(i,j) = f(idx);

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

fprintf('Pass rates: A %.0f%%   B %.0f%%   C %.0f%%   ALL %.1f%%\n', ...
        100*mean(passA(:)), 100*mean(passB(:)), ...
        100*mean(passC(:)), 100*mean(feasible(:)));

if ~any(feasible(:))
    fprintf(['\nNOTHING PASSES. Read the individual rates above before\n' ...
             'touching any criterion - they tell you which one is\n' ...
             'binding. If C is the culprit, suspect the criterion\n' ...
             '(FLAT_TOL_DB is an assumption, not a published spec)\n' ...
             'before suspecting the model.\n']);
end

%% ------------------------------------------------------------------
%  Plots
%  ------------------------------------------------------------------
[RR, LL] = meshgrid(R_vec, L1_vec*1e3);

figure('Name','Velocity envelope','Position',[100 100 560 380]);
loglog(f, v_env, 'LineWidth', 2); grid on; hold on;
xline(f_corner_lo, 'k--'); xline(f_corner_hi, 'k--');
xlabel('Frequency (Hz)'); ylabel('Peak velocity (cm/s)');
title('Recorded velocity envelope (three physical limits)');

figure('Name','Cartridge parameter feasibility','Position',[100 100 1100 800]);

subplot(2,3,1);
contourf(RR, LL, f_hf/1e3, 20, 'LineColor','none'); hold on;
contour(RR, LL, f_hf, [F_MIN_HF F_MIN_HF], 'w-', 'LineWidth', 2);
colorbar; xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title('A: HF resonance (kHz)');

subplot(2,3,2);
contourf(RR, LL, log10(peak_force), 20, 'LineColor','none'); hold on;
contour(RR, LL, peak_force, [VTF_MAX VTF_MAX], 'w-', 'LineWidth', 2);
cb = colorbar; cb.Label.String = 'log_{10}(gram-force)';
xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title(sprintf('B: peak force demand, limit %.1f gf', VTF_MAX));

subplot(2,3,3);
contourf(RR, LL, ripple_dB, 20, 'LineColor','none'); hold on;
contour(RR, LL, ripple_dB, [FLAT_TOL_DB FLAT_TOL_DB], 'w-', 'LineWidth', 2);
colorbar; xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title('C: response ripple (dB p-p)');

% Margin plots - how much each criterion passes/fails by. Far more
% informative than a binary map when nothing passes.
subplot(2,3,4);
margB = 20*log10(VTF_MAX ./ peak_force);
contourf(RR, LL, margB, 20, 'LineColor','none'); hold on;
contour(RR, LL, margB, [0 0], 'w-', 'LineWidth', 2);
colorbar; xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title('B margin (dB, >0 passes)');

subplot(2,3,5);
contourf(RR, LL, f_at_peak/1e3, 20, 'LineColor','none');
colorbar; xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title('Frequency of peak force demand (kHz)');

subplot(2,3,6);
imagesc(R_vec, L1_vec*1e3, double(feasible));
set(gca,'YDir','normal'); clim([0 1]);
colormap(gca, [0.88 0.88 0.88; 0.15 0.55 0.35]);
xlabel('R (dyne-s/cm)'); ylabel('L1 (mg)');
title('Feasible (all three)');

%% ------------------------------------------------------------------
%  Report
%  ------------------------------------------------------------------
if any(feasible(:))
    [iF, jF] = find(feasible);
    fprintf('\nFeasible region spans:\n');
    fprintf('  L1: %.3f to %.3f mg\n', min(L1_vec(iF))*1e3, max(L1_vec(iF))*1e3);
    fprintf('  R : %.1f to %.1f dyne-s/cm\n', min(R_vec(jF)), max(R_vec(jF)));
    fprintf('\nCentroid (convenience only - meaningless if the region\n');
    fprintf('is L-shaped or disconnected, so look at the map first):\n');
    fprintf('  L1 = %.3f mg,  R = %.1f\n', ...
            mean(L1_vec(iF))*1e3, mean(R_vec(jF)));
end

%% ------------------------------------------------------------------
%  ZENER BEARING - clarifying the V1 comment, which was wrong
%  ------------------------------------------------------------------
% V1's comment implied you were waiting on a datasheet value you did not
% have. You are not, and C2 was never a swept parameter.
%
% C2 above IS the datasheet number - the DYNAMIC compliance at 100 Hz.
% The sweep only ever estimates L1 and R.
%
% The datasheet ALSO gives a STATIC compliance, and the simple C2+R pair
% has nowhere to put it. The Zener (standard linear solid) model does,
% and stays rational so the whole tf/lsim/series toolchain still works:
%
%     Ceff(s) = Cinf + (C0 - Cinf)/(1 + s*tau)
%     Zbear(s) = 1/(s*Ceff(s))
%
%   C0   <- datasheet STATIC compliance
%   Cinf <- datasheet DYNAMIC compliance (100 Hz, ~ the HF asymptote)
%   tau  <- the only genuinely free parameter
%
% This REPLACES C2 and R together - relaxation IS the damping, so
% keeping both double-counts one physical mechanism. The sweep then
% becomes (L1, tau) rather than (L1, R), with better physical grounding
% and one fewer unconstrained parameter.
%
% Sanity check on the current model: it already has implicit frequency-
% dependent compliance, Ceff(s) = C2/(1 + s*R*C2), cornering at
% 1/(R*C2) = 159 Hz for R = 40. Suspiciously close to AT's 100 Hz
% measurement point, which suggests the article's values are physically
% sane. What the two-parameter form gets wrong is driving compliance to
% zero at HF rather than to a finite glassy plateau - exactly what the
% third parameter fixes.