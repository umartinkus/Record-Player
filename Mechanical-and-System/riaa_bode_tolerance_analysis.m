%RIAA_BODE_TOLERANCE_ANALYSIS
%   Bode plot of the RIAA amplifier (see riaa_amp_tf.m) with an
%   uncertainty band computed by first-order error propagation through
%   the transfer function, given assumed per-component tolerances.
%   Also overlays the ideal IEC-amended RIAA target curve for
%   reference, so you can see both "how far is the nominal design from
%   ideal" and "how much does part tolerance spread that."
%
%   METHOD (error propagation, not Monte Carlo):
%   At each frequency, the magnitude response in dB is a function of
%   six component values, dB = f(R0,C0,R1,C1,R2,C2). For each
%   component we estimate the local sensitivity dDB/dp_i numerically
%   (central finite difference - perturb the part up and down by a
%   tiny fraction and see how much the dB response moves). Given each
%   component's assumed tolerance (treated here as a 1-sigma spread),
%   the standard first-order uncertainty-propagation formula combines
%   them:
%
%       sigma_dB(f) = sqrt( sum_i ( dDB/dp_i * sigma_i )^2 )   [RSS]
%
%   This assumes the component errors are independent and small enough
%   that the linear (first-order Taylor) approximation holds - true
%   here since tolerances are a few percent, not tens of percent.
%   RSS is the "realistic" band: how far a typical real board is
%   likely to land, since it's very unlikely every component errs in
%   the same direction simultaneously. We also compute the "worst
%   case" bound (every component simultaneously at its tolerance limit,
%   in whichever direction hurts most) as a linear sum instead of RSS -
%   useful to know as a paranoid upper bound, but expect real boards to
%   sit well inside it.
%
%   Requires: riaa_amp_tf.m on the MATLAB path, Control System Toolbox.

clear; clc; close all;

%% Nominal component values - AS-BUILT standard-value parts, not the
%  book's exact numbers (riaa_amp_tf's zero-argument default uses the
%  book's exact values; here we use the nearest standard E96/E24 parts
%  actually going on the board, per the tolerance/sourcing discussion).
% Reference: Douglas Self's Small Signal Audio Design (2024)
R0 = 1e3;      C0 = 7.96e-6;    % electrolytic, loose tolerance regardless of nominal
R1 = 511e3;    C1 = 6.2e-9;
R2 = 43.2e3;   C2 = 1.8e-9;     % update this (and its tolerance below) if you
                                % combine two parts for tighter accuracy here

p0    = [R0, C0, R1, C1, R2, C2];
names = {'R0','C0','R1','C1','R2','C2'};

%% Per-component tolerance (fractional, e.g. 0.01 = 1%), treated as 1-sigma
% 1% metal-film resistors, 5% precision film/C0G caps for the RIAA
% timing network, 20% for the electrolytic C0 (typical for electrolytics
% - this breakpoint is the subsonic pole, not one of the core RIAA
% timing constants, so it's far less sensitive to begin with).
tol = [0.01, 0.20, 0.01, 0.05, 0.01, 0.05];   % [R0 C0 R1 C1 R2 C2]

%% Frequency sweep
f = logspace(log10(10), log10(100e3), 400);   % 10 Hz to 100 kHz
w = 2*pi*f;                                    % rad/s, what bode() wants

%% Nominal response
G0 = riaa_amp_tf(p0(1), p0(2), p0(3), p0(4), p0(5), p0(6));
[mag0, phase0] = bode(G0, w);
mag0_dB   = squeeze(20*log10(mag0));
phase0deg = squeeze(phase0);

%% Sensitivity via central-difference finite differences
n = numel(p0);
dMagdP   = zeros(numel(f), n);   % d(dB)/d(param), one column per component
dPhasedP = zeros(numel(f), n);   % d(deg)/d(param)

frac_step = 1e-4;  % 0.01% nudge - small enough for a clean numerical
                    % derivative, unrelated to the component's real tolerance

for i = 1:n
    pPlus  = p0;  pPlus(i)  = p0(i) * (1 + frac_step);
    pMinus = p0;  pMinus(i) = p0(i) * (1 - frac_step);
    h = pPlus(i) - pMinus(i);

    Gp = riaa_amp_tf(pPlus(1),  pPlus(2),  pPlus(3),  pPlus(4),  pPlus(5),  pPlus(6));
    Gm = riaa_amp_tf(pMinus(1), pMinus(2), pMinus(3), pMinus(4), pMinus(5), pMinus(6));

    [magP, phaseP] = bode(Gp, w);
    [magM, phaseM] = bode(Gm, w);

    dMagdP(:,i)   = (squeeze(20*log10(magP)) - squeeze(20*log10(magM))) / h;
    dPhasedP(:,i) = (squeeze(phaseP)         - squeeze(phaseM))         / h;
end

%% Propagate uncertainty
sigma = tol .* p0;   % 1-sigma of each component, in its own native units (ohms/farads)

% RSS (statistical / "realistic spread across real boards") combination.
% Written as a matrix-vector product: for each frequency row, this sums
% (dMagdP_i)^2 * sigma_i^2 over the 6 components in one shot.
magVar        = (dMagdP  .^ 2) * (sigma' .^ 2);
sigma_mag_dB  = sqrt(magVar);
phaseVar      = (dPhasedP.^ 2) * (sigma' .^ 2);
sigma_phase_deg = sqrt(phaseVar);

% Worst-case (linear sum / "every component simultaneously unlucky") bound.
worst_mag_dB   = abs(dMagdP)   * sigma';
worst_phase_deg = abs(dPhasedP) * sigma';

%% Ideal IEC-amended RIAA target curve, for comparison only
T1 = 3180e-6; T2 = 318e-6; T3 = 75e-6; T4 = 7950e-6;   % 50Hz, 500Hz, 2122Hz, 20Hz(IEC)
s_jw  = 1i*w;
Hriaa = (1 + s_jw*T2) ./ ((1 + s_jw*T1) .* (1 + s_jw*T3) .* (1 + s_jw*T4));
Hriaa_dB = 20*log10(abs(Hriaa));

% Normalize both curves to 0 dB at 1 kHz for a fair shape comparison
[~, idx1k]      = min(abs(f - 1000));
Hriaa_dB_norm   = Hriaa_dB - Hriaa_dB(idx1k);
mag0_dB_norm    = mag0_dB  - mag0_dB(idx1k);

%% --- Plot: magnitude Bode with uncertainty band ---
figure('Name', 'RIAA Amplifier - Bode with Tolerance Error Propagation', ...
       'Position', [100 100 900 700]);

subplot(2,1,1);
k_sigma = 3;   % show +/-3 sigma (~99.7% coverage if errors are ~Gaussian)
upper3 = mag0_dB_norm + k_sigma*sigma_mag_dB;
lower3 = mag0_dB_norm - k_sigma*sigma_mag_dB;

fill([f, fliplr(f)], [upper3', fliplr(lower3')], [0.85 0.90 1.0], ...
     'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
     'DisplayName', sprintf('\\pm%d\\sigma (RSS, realistic spread)', k_sigma));
hold on;
plot(f, mag0_dB_norm + worst_mag_dB, 'r--', 'DisplayName', 'Worst case (all parts unlucky)');
plot(f, mag0_dB_norm - worst_mag_dB, 'r--', 'HandleVisibility', 'off');
plot(f, Hriaa_dB_norm, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Ideal IEC RIAA target');
plot(f, mag0_dB_norm, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Nominal (as-built values)');

% Literal error bars at standard reference frequencies
f_ref   = [20 50 100 500 1000 2122 10000 20000];
mag_ref = interp1(f, mag0_dB_norm,   f_ref);
sig_ref = interp1(f, sigma_mag_dB,   f_ref);
errorbar(f_ref, mag_ref, k_sigma*sig_ref, 'o', 'Color', [0 0.45 0], ...
    'MarkerFaceColor', [0 0.45 0], 'LineWidth', 1.2, ...
    'DisplayName', sprintf('\\pm%d\\sigma at reference points', k_sigma));

set(gca, 'XScale', 'log');
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB, normalized to 0 dB @ 1 kHz)');
title('RIAA Amplifier Frequency Response with Component Tolerance Uncertainty');
legend('Location', 'southwest');
xlim([10 100e3]);

subplot(2,1,2);
upperP = phase0deg + k_sigma*sigma_phase_deg;
lowerP = phase0deg - k_sigma*sigma_phase_deg;
fill([f, fliplr(f)], [upperP', fliplr(lowerP')], [0.85 0.90 1.0], ...
     'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
     'DisplayName', sprintf('\\pm%d\\sigma (RSS)', k_sigma));
hold on;
plot(f, phase0deg, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Nominal');
set(gca, 'XScale', 'log');
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
title('Phase Response with Tolerance Uncertainty');
legend('Location', 'best');
xlim([10 100e3]);

%% Report which component dominates the error, by frequency band
fprintf('\nDominant error contributor by frequency band (largest avg |dMag/dp * sigma|):\n');
bands = [20 100; 100 1000; 1000 10000; 10000 20000];
for b = 1:size(bands,1)
    idxBand = f >= bands(b,1) & f <= bands(b,2);
    contrib = mean(abs(dMagdP(idxBand,:)) .* sigma, 1);
    [maxContrib, worstIdx] = max(contrib);
    fprintf('  %6.0f - %6.0f Hz: %s  (avg contribution %.4f dB)\n', ...
        bands(b,1), bands(b,2), names{worstIdx}, maxContrib);
end
