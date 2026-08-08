L1 = 0.80e-3; C2 = 25e-6; C1 = 0.0506e-6; R = 40;
[Hvel, Hforce] = cartridge_model_mech(L1, C1, C2, R);

fprintf('Series (force notch):  %.0f Hz\n', 1/(2*pi*sqrt(L1*C2)));
fprintf('Parallel (force peak): %.0f Hz\n', 1/(2*pi*sqrt(L1*C1)));

figure; bode(Hvel, Hforce, {2*pi*20, 2*pi*50e3});
grid on;
legend('Velocity response','Tracking force');