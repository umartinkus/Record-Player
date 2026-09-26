"""Apply a RIAA playback filter to a recorded stereo signal.

Python port of DSP.m.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.io import wavfile


def write_wav16(path, x, fs):
    """Mimic MATLAB audiowrite(..., 'BitsPerSample', 16) for float data:
    samples are clipped to [-1, 1] and scaled to int16."""
    x = np.clip(x, -1.0, 1.0)
    wavfile.write(path, fs, np.round(x * 32767).astype(np.int16))


# Import left and right channel of audio signal
y = np.loadtxt('RecData.csv', delimiter=',')
yNorm = y / np.max(np.abs(y), axis=0)  # normalize original sig
Fs = 80000
t = np.arange(1, len(y) + 1) / Fs
A0 = 0.7e4  # scaling factor to place 0dB point at 1kHz

# Define RIAA Playback transfer function
RIAA_Filter = signal.TransferFunction(A0 * np.array([1, 500, 0]),
                                      [1, 2070, 141000, 2000000])
f_bode = np.logspace(0, 5, 1000)  # Hz
_, h = signal.freqs(RIAA_Filter.num, RIAA_Filter.den, worN=2 * np.pi * f_bode)

bodeFig, (ax_mag, ax_ph) = plt.subplots(2, 1, sharex=True)
ax_mag.semilogx(f_bode, 20 * np.log10(np.abs(h)))
ax_mag.set_ylabel('Magnitude (dB)')
ax_mag.set_title('Bode Plot of RIAA Transfer Function')
ax_mag.grid(True, which='both')
ax_ph.semilogx(f_bode, np.degrees(np.unwrap(np.angle(h))))
ax_ph.set_ylabel('Phase (deg)')
ax_ph.set_xlabel('Frequency (Hz)')
ax_ph.grid(True, which='both')

# Send input data through RIAA Filter.
# Equivalent to MATLAB lsim (first-order hold on the input): discretize
# with FOH at the sample rate, then run a fast IIR filter.
bd, ad, _ = signal.cont2discrete((RIAA_Filter.num, RIAA_Filter.den),
                                 1 / Fs, method='foh')
bd = np.squeeze(bd)
yFL = signal.lfilter(bd, ad, y[:, 0])
yFL = yFL / np.max(np.abs(yFL))  # normalize left channel
yFR = signal.lfilter(bd, ad, y[:, 1])
yFR = yFR / np.max(np.abs(yFR))  # normalize right channel

# Write data to a file
write_wav16('RawOutput.wav', yNorm, Fs)
write_wav16('TestOutput.wav', y, Fs)

# Plot PSD of the raw and filtered signal
# (MATLAB pwelch defaults for an integer window: Hamming, 50% overlap)
filterFig = plt.figure()
N = 512
F, Px = signal.welch(np.column_stack([yFL, yFR, y]), fs=Fs, window='hamming',
                     nperseg=N, noverlap=N // 2, nfft=N, axis=0)
plt.plot(np.log10(F), 20 * np.log10(Px))  # Plots the power spectrum
# scaling F by 1000 will represent frequency in kHz
plt.xlabel('Frequency (log10(Hz)) ')
plt.ylabel('Power Spectral Density (in dB) ')
plt.legend(['Filtered LCH', 'Filtered RCH', 'Raw LCH', 'Raw RCH'])

plt.show()
