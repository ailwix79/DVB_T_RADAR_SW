% Coarse Symbol Synchronization With Rational Delay
% First OFDM Symbol Position Detection
% Montecarlo Simulations of estimation algorithm

% Author: Alejandro Manuel López Gómez (simulation code)
% Author: Carlos García de la Cueva (CSS algorithm and dame_dvbt_bb_signal)

% Parameters
% - Signal parameters. Used to create the DVB-T signal. With the current
% numbers only one OFDM symbol is created, to which later a delay is
% introduced.
% - Function parameters. These define minimum and maximum SNR values (used for the noise introduction), the
% desired number of simulations for each SNR value, the signal
% smoothing / averaging values, maximum delay allowed for the signal and
% the sample error threshold.

% TODO. Explore other signal smoothing methods (peak envelope, different
% filter types (Savitzky-Golay)...)

% CHANGES TO PREVIOUS VERSION
% - Introduced the dame_dvbt_bb_signal function for singla generation
% - Changed the noise generator to awgn Matlab function
% - Added the possibility of plotting more than one averaged signal for
% conclusion extraction purposes

clear all;
close all;

% Signal Parameters
BW = 8;                                 % Channel selection 5,6,7,8 (Mhz)
fs = (64/7)*1e6;
tx_mode = '8K';                         % Transmitter mode '2K','4K','8K', different number of carriers and Tu values for each mode...
frame_offset = 0;                    % Frame Offset within a superframe (PREGUNTAR)
guard = 1/32;                           % Guard interval length (fraction of T, page 33 table 14)
T_symb = 8192;
CP = guard*T_symb;
mod_type = '64-QAM';                    % Symbol modulation 'QPSK','16-QAM','64-QAM'
alpha = 2;                              % Normalization factor for the modulation of the OFDM symbol (page 27)

% Function Parameters
min_SNR = -30;                          % Minimum SNR value for simulations
max_SNR = 30;                           % Maximum SNR value
n_simul = 1000;                         % Number of simulations per SNR value
avg_value_set = [1,4,8,16,32];        % Signal Averaging values (no size limitation)
max_delay = 32;                        % Maximum delay possible
max_accepted_rmse = 5;                  % Maximum RMSE sample error accepted

snr_set = min_SNR:1:max_SNR;                                                    % SNR value range
E = zeros(length(avg_value_set)*length(snr_set),n_simul);                       % Non averaged signal error

tic

for i=1:length(avg_value_set)
    
    T = ((avg_value_set(i)+1)*(T_symb+CP))/fs;
    
    for ii=1:length(snr_set)
        
        % DVB-T Signal Generation (One new signal per SNR value)
        [s_tx,~,~,~] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);
        
        for iii=1:n_simul
                if floor(100*iii/(n_simul)) > floor(100*(iii-1)/(n_simul))
                    clc;
                    fprintf('Realizando simulaciones para SNR: %d dB... (%3.0f %%)\n',snr_set(ii),floor(100*iii/(n_simul)));
                    fprintf('Valor de promediado: %d (%d / %d)\n',avg_value_set(i),i,length(avg_value_set));
                end
                % Add AWGN
                s_tx_w_n = awgn(s_tx,snr_set(ii),'measured')'; 

                % Introduce delay
                n_tram = max_delay*rand;
                lfft = length(s_tx_w_n) + max_delay - 1;
                s_tx_w_n = ifft(fft(s_tx_w_n,lfft).*exp(-1j*2*pi*n_tram*(1/(lfft)).*(0:lfft-1)),lfft);
                s_tx_w_n = s_tx_w_n';

                % Coarse Symbol Estimation Algorithm
                s_x = conj([s_tx_w_n; zeros(T_symb,1)]).*[zeros(T_symb,1); s_tx_w_n];       % Correlation function
                s_x = conv(s_x, ones(CP,1));                                                % Signal smoothing
                s_x = s_x(1:avg_value_set(i)*(T_symb+CP));                                  % Adjust signal length to be a multiple of the averaging value
                A = reshape(s_x,T_symb+CP,length(s_x)/(T_symb+CP));
                estimate = sum(abs(A).^2,2);                               

                % Estimation extraction
                [~,peak_index] = max(estimate);

                % Error calculation
                if (peak_index > (T_symb+CP)/2)
                    peak_index = peak_index - (T_symb+CP);
                end

                E(ii+(i-1)*length(snr_set),iii) = n_tram - peak_index;
        end
    end
end

% RMS calculation (column dimension due to how data has been allocated)
for i=1:length(avg_value_set)
    RMS_E(((i-1)*length(snr_set) + 1):i*length(snr_set),:) = rms(E(((i-1)*length(snr_set) + 1):i*length(snr_set),:),2);
end

% RMS error plot
figure;
hold on;
grid on;

for i=1:length(avg_value_set)
    plot(snr_set,RMS_E(((i-1)*length(snr_set) + 1):i*length(snr_set),:));
end
    
yline(max_accepted_rmse,'--','Threshold','LineWidth',2);
legend(append('Averaging ',cellstr(int2str(avg_value_set'))));
title(['RMS Error (rational delay). Number of simulations per SNR value = ',num2str(n_simul)]);
ylabel('RMS Sample Error');
xlabel('SNR (dB)');
set(gca,'Yscale','log');

toc
