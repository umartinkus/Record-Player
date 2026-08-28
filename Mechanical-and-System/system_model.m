% L1_mH = 0.80e-3; C2_uF = 25e-6; C1_uF = 0.0506e-6; R_ohm = 40;
% L1_mH = 0.8e-3; C1_uF = 0.0506e-6; C2_uF = 17e-6; R_ohm = 40;
% [Hvel, Hforce] = cartridge_tf(L1_mH, C1_uF, C2_uF, R_ohm);
% 
% fprintf('Series (force notch):  %.0f Hz\n', 1/(2*pi*sqrt(L1_mH*C2_uF)));
% fprintf('Parallel (force peak): %.0f Hz\n', 1/(2*pi*sqrt(L1_mH*C1_uF)));
% 
% figure; bode(Hvel, Hforce, {2*pi*20, 2*pi*50e3});
% % figure; pzplot(Hvel, Hforce);
% grid on;
% legend('Velocity response','Tracking force');

sys = cartridge_tf_mech*cartridge_tf_elec;
bode(sys)