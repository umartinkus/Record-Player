% Import left and right channel of audio signal
y = csvread("RecData.csv");
yNorm = [y(:,1)/max(abs(y(:,1))) y(:,2)/max(abs(y(:,2)))]; % normalize original sig
Fs = 80000;
t = (1:length(y))/Fs;
A0 = 0.7e4; % scaling factor to place 0dB point at 1kHz

% Define RIAA Playback transfer function
RIAA_Filter = A0*tf([1 500 0],[1 2070 141000 2000000]);
P = bodeoptions;
P.FreqUnits = 'Hz';
bodeplot(RIAA_Filter,P);
bodeFig = gca;
title(bodeFig,'Bode Plot of RIAA Transfer Function')


% Send input data through RIAA Filter
yFL = lsim(RIAA_Filter,y(:,1),t);
yFL = yFL/max(abs(yFL)); % normalize left channel
yFR = lsim(RIAA_Filter,y(:,2),t);
yFR = yFR/max(abs(yFR)); %normalize right channel


% Write data to a file
audiowrite("RawOutput.wav",yNorm,Fs,'BitsPerSample',16)
audiowrite("TestOutput.wav",y,Fs,'BitsPerSample',16)

% Plot PSD of the raw and filtered signal
filterFig = figure;
N = 512; 
[Px, F]=pwelch([yFL yFR y], N, [], N, Fs); 
plot(log10(F), 20*log10(Px)); %Plots the power spectrum 
%scaling F by 1000 will represent frequency in kHz 
xlabel('Frequency (log10(Hz)) '); 
ylabel('Power Spectral Density (in dB) '); 
legend('Filtered LCH','Filtered RCH','Raw LCH','Raw RCH')
 