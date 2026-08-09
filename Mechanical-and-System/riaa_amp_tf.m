function G = riaa_amp_tf(R0, C0, R1, C1, R2, C2)
%RIAA_AMP_TF  Laplace-domain transfer function of the single-stage,
%   non-inverting active RIAA amplifier (Douglas Self, "Small Signal
%   Audio Design"). This models the OP-AMP GAIN STAGE ONLY: input at
%   the op-amp's non-inverting (+) terminal, output at the op-amp
%   output. It does NOT include the cartridge's source impedance,
%   47k/pF loading network, or cable capacitance - those belong in a
%   separate front-end block when you assemble the full system model.
%
%   G = RIAA_AMP_TF() returns the transfer function using the exact
%       reference component values from Self's book:
%           R0 = 1e3        C0 = 7.96e-6
%           R1 = 513.5e3    C1 = 6.193e-9
%           R2 = 42.86e3    C2 = 1.75e-9
%
%   G = RIAA_AMP_TF(R0,C0,R1,C1,R2,C2) uses custom component values
%       (SI units throughout: ohms, farads). Call this with perturbed
%       values for tolerance sweeps, Monte Carlo analysis, or to model
%       the as-built board using real (standard-value) parts instead
%       of the book's exact numbers.
%
%   Circuit (non-inverting op-amp topology):
%
%       Gain(s) = 1 + Zf(s) / Zg(s)
%
%       Zg(s) = R0 + 1/(s*C0)                      ground leg, from the
%                                                    inverting (-) input
%       Zf(s) = R1/(1+s*R1*C1) + R2/(1+s*R2*C2)     feedback network,
%                                                    output to (-) input
%
%   R0*C0 sets the subsonic/DC-gain-unity time constant (~20 Hz, the
%   IEC amendment pole - it's what keeps DC gain at unity and rolls off
%   turntable rumble, built into the network rather than a separate
%   stage). R1*C1 = 3180us (50 Hz pole). R2*C2 = 75us (2122 Hz pole).
%   The 500 Hz breakpoint is not a direct RC product - it emerges from
%   the interaction of the two feedback branches (see Lipshitz, "On
%   RIAA Equalization Networks", JAES 1979, for the full derivation).
%
%   Returns G as a Control System Toolbox `tf` object (Laplace/s
%   domain). Requires the Control System Toolbox (bode, freqresp,
%   evalfr, minreal, etc. all work directly on the returned object).
%
%   Example:
%       G = riaa_amp_tf();
%       bode(G); grid on;
%
%       % gain at 1 kHz, in dB
%       g1k_dB = 20*log10(abs(evalfr(G, 1i*2*pi*1000)));
%
%       % swap in as-built standard-value parts instead of the book's
%       % exact values
%       G_asbuilt = riaa_amp_tf(1e3, 8.2e-6, 511e3, 6.2e-9, 43.2e3, 1.8e-9);
%
%   See also: riaa_bode_tolerance_analysis.m (separate analysis script
%   that calls this function to build a Bode plot with a component-
%   tolerance uncertainty band via error propagation).

    if nargin == 0
        R0 = 1e3;      C0 = 7.96e-6;
        R1 = 513.5e3;  C1 = 6.193e-9;
        R2 = 42.86e3;  C2 = 1.75e-9;
    elseif nargin ~= 6
        error('riaa_amp_tf:nargin', ...
            ['Call with zero arguments (book defaults) or all six: ', ...
             'riaa_amp_tf(R0,C0,R1,C1,R2,C2).']);
    end

    s = tf('s');

    Zg = R0 + 1/(s*C0);
    Z1 = R1 / (1 + s*R1*C1);
    Z2 = R2 / (1 + s*R2*C2);
    Zf = Z1 + Z2;

    G = 1 + Zf/Zg;
    G = minreal(G);   % clean up any near-common pole/zero pairs from the tf algebra
end
