% Caracterización de la función de ambiguedad

%% PARAMETROS INICIALES

clc;
clear all;

BW = 8;                                 % Channel selection 5,6,7,8 (Mhz)
n_symb = 128;                           % Number of symbols
fs = (64/7)*1e6;                        % Sampling frequency (Hz)
tx_mode = '8K';                         % Transmitter mode '2K','4K','8K', different number of carriers and Tu values for each mode...
frame_offset = 0;                       % Frame offset
guard = 1/32;                           % Guard interval length (fraction of T, page 33 table 14)
T_symb = 8192;                          % Useful symbol time in samples
CP = guard*T_symb;                      % Guard length in samples
mod_type = '64-QAM';                    % Symbol modulation 'QPSK','16-QAM','64-QAM'
alpha = 2;                              % Normalization factor for the modulation of the OFDM symbol (page 27)
T = (n_symb*(T_symb + CP))/fs;          % Signal duration (in seconds)
N = T*fs;                               % Signal length (it will serve as the CPI in this case)
PRF = 1/T;


%% NO ENVENTANADO

[s_x,~] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);

[afmagdel,delay] = ambgfun(s_x,fs,PRF,"Cut","Doppler");
[afmagdop,doppler] = ambgfun(s_x,fs,PRF,"Cut","Delay");

figure;
sgtitle("DVB-T Ambiguity function analysis")
subplot(2,2,1)
plot(delay,20*log10(afmagdel),'r')
grid on;
title("Delay")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,2)
plot(doppler,20*log10(afmagdop),'b')
grid on;
title("Doppler shift")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")

subplot(2,2,3)
plot(delay,20*log10(afmagdel),'r')
xlim([-1e-5,1e-5]);
grid on;
title("Delay (Zero delay closeup)")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,4)
plot(doppler,20*log10(afmagdop),'b')
xlim([-100,100]);
grid on;
title("Doppler shift (Zero Doppler closeup)")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")

%% ENVENTANADO EN EL TIEMPO

[s_x,~] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);

s_x = s_x .* chebwin(length(s_x));

[afmagdel,delay] = ambgfun(s_x,fs,PRF,"Cut","Doppler");
[afmagdop,doppler] = ambgfun(s_x,fs,PRF,"Cut","Delay");

figure;
sgtitle("DVB-T Ambiguity function analysis (TIME windowed)")
subplot(2,2,1)
plot(delay,20*log10(afmagdel),'r')
grid on;
title("Delay")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,2)
plot(doppler,20*log10(afmagdop),'b')
grid on;
title("Doppler shift")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")

subplot(2,2,3)
plot(delay,20*log10(afmagdel),'r')
xlim([-1e-5,1e-5]);
grid on;
title("Delay (Zero delay closeup)")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,4)
plot(doppler,20*log10(afmagdop),'b')
xlim([-100,100]);
grid on;
title("Doppler shift (Zero Doppler closeup)")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")


%% ENVENTANADO EN FRECUENCIA

[s_x,~] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);

s_x = ifft(ifftshift(fftshift(fft(s_x,length(s_x))).* chebwin(length(s_x))),length(s_x));

[afmagdel,delay] = ambgfun(s_x,fs,PRF,"Cut","Doppler");
[afmagdop,doppler] = ambgfun(s_x,fs,PRF,"Cut","Delay");

figure;
sgtitle("DVB-T Ambiguity function analysis (FREQUENCY windowed)")
subplot(2,2,1)
plot(delay,20*log10(afmagdel),'r')
grid on;
title("Delay")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,2)
plot(doppler,20*log10(afmagdop),'b')
grid on;
title("Doppler shift")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")

subplot(2,2,3)
plot(delay,20*log10(afmagdel),'r')
xlim([-1e-5,1e-5]);
grid on;
title("Delay (Zero delay closeup)")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,4)
plot(doppler,20*log10(afmagdop),'b')
xlim([-100,100]);
grid on;
title("Doppler shift (Zero Doppler closeup)")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")


%% INTERPOLACION Y ENVENTANADO EN EL TIEMPO

[s_x,~] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);

L = 2;
s_x = interp(s_x,L);
s_x = s_x .* L .* [zeros((length(s_x)-length(s_x)/L)/2,1) ; ones(length(s_x)/L,1) ; zeros((length(s_x)-length(s_x)/L)/2,1)];
s_x(s_x == 0) = [];

s_x = s_x .* chebwin(length(s_x));

[afmagdel,delay] = ambgfun(s_x,fs,PRF,"Cut","Doppler");
[afmagdop,doppler] = ambgfun(s_x,fs,PRF,"Cut","Delay");

figure;
sgtitle("DVB-T Ambiguity function analysis (TIME win and L = 2)")
subplot(2,2,1)
plot(delay,20*log10(afmagdel),'r')
grid on;
title("Delay")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,2)
plot(doppler,20*log10(afmagdop),'b')
grid on;
title("Doppler shift")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")

subplot(2,2,3)
plot(delay,20*log10(afmagdel),'r')
xlim([-1e-5,1e-5]);
grid on;
title("Delay (Zero delay closeup)")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,4)
plot(doppler,20*log10(afmagdop),'b')
xlim([-100,100]);
grid on;
title("Doppler shift (Zero Doppler closeup)")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")

%% INTERPOLACIÓN Y ENVENTANADO EN FRECUENCIA

[s_x,~] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);

L = 2;
s_x = interp(s_x,L);
s_x = s_x .* L .* [zeros((length(s_x)-length(s_x)/L)/2,1) ; ones(length(s_x)/L,1) ; zeros((length(s_x)-length(s_x)/L)/2,1)];
s_x(s_x == 0) = [];

s_x = ifft(ifftshift(fftshift(fft(s_x,length(s_x))).* chebwin(length(s_x))),length(s_x));

[afmagdel,delay] = ambgfun(s_x,fs,PRF,"Cut","Doppler");
[afmagdop,doppler] = ambgfun(s_x,fs,PRF,"Cut","Delay");

figure;
sgtitle("DVB-T Ambiguity function analysis (FREQUENCY win and L = 2)")
subplot(2,2,1)
plot(delay,20*log10(afmagdel),'r')
grid on;
title("Delay")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,2)
plot(doppler,20*log10(afmagdop),'b')
grid on;
title("Doppler shift")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")

subplot(2,2,3)
plot(delay,20*log10(afmagdel),'r')
xlim([-1e-5,1e-5]);
grid on;
title("Delay (Zero delay closeup)")
xlabel("Delay (seconds)")
ylabel("Normalised Power level (dB)")

subplot(2,2,4)
plot(doppler,20*log10(afmagdop),'b')
xlim([-100,100]);
grid on;
title("Doppler shift (Zero Doppler closeup)")
xlabel("Doppler frequency (Hz)")
ylabel("Normalised Power level (dB)")

