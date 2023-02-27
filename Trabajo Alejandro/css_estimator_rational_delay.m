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
% TODO. Introduce a more complex DVB_T signal (not just only one symbol)

% CHANGES TO PREVIOUS VERSIONS
% - Introduced the dame_dvbt_bb_signal function for singla generation
% - Changed the noise generator to awgn Matlab function
% - Added the possibility of plotting more than one averaged signal for
% conclusion extraction purposes

clear all;
close all;

% Signal Parameters
BW = 5;                                 % Channel selection 5,6,7,8 (Mhz)
T = 1e-5;                               % Symbol duration
tx_mode = '2K';                         % Transmitter mode '2K','4K','8K', different number of carriers and Tu values for each mode...
frame_offset = T/16;                    % Frame Offset within a superframe (PREGUNTAR)
guard = 1/32;                           % Guard interval length (fraction of T, page 33 table 14)
mod_type = '64-QAM';                    % Symbol modulation 'QPSK','16-QAM','64-QAM'
alpha = 2;                              % Normalization factor for the modulation of the OFDM symbol (page 27)

% Function Parameters
min_SNR = -30;                          % Minimum SNR value for simulations
max_SNR = 30;                           % Maximum SNR value
n_simul = 1000;                         % Number of simulations per SNR value
avg_value_set = [8,16,32,34,36,38,40];  % Signal Averaging values (no size limitation)
max_delay = 512;                        % Maximum delay possible
max_accepted_rmse = 4;                  % Maximum RMSE sample error accepted

snr_set = min_SNR:1:max_SNR;                                                    % SNR value range
E = zeros(length(avg_value_set)*length(snr_set) + length(snr_set),n_simul);     % Non averaged signal error

tic
for i=1:length(snr_set)
        
        % State Check
        if floor(100*i/length(snr_set)) > floor(100*(i-1)/(length(snr_set)))
          clc;
          fprintf('In process ... (%3.0f %%)\n',floor(100*i/(length(snr_set))));
        end
        
        % DVB-T Signal Generation (One new signal per SNR value)
        [s_tx,T_symb] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);
        CP = length(s_tx)-T_symb;
    
        for ii=1:n_simul
        
            % Add AWGN
            s_tx_w_n = awgn(s_tx,snr_set(i),'measured')'; 
%             test = s_tx_w_n;
            
            % Introduce delay
            n_tram = max_delay*rand;
            lfft = length(s_tx_w_n) + max_delay - 1;
            s_tx_w_n = ifft(fft(s_tx_w_n,lfft).*exp(-1j*2*pi*n_tram*(1/(lfft)).*(0:lfft-1)),lfft);
            s_tx_w_n = s_tx_w_n';
            
            % Verify delay
%             figure;
%             hold on;
%             plot(abs(test));
%             plot(abs(s_tx_w_n));
            
            % Check delay (WIP)
%             d1 = finddelay(test,s_tx_w_n);
            
            % Coarse Symbol Estimation Algorithm
            for iii=1:length(avg_value_set)
                s_x = conj([s_tx_w_n; zeros(T_symb,1)]).*[zeros(T_symb,1); s_tx_w_n];       % Correlation function
                s_x = conv(s_x, ones(CP,1));                                                % Signal smoothing

                s_x = [s_x ; zeros(T_symb+CP-mod(length(s_x),T_symb+CP),1)];                % Adjust signal length for FFT.
                A = reshape(s_x,T_symb+CP,length(s_x)/(T_symb+CP));
                estimate = sum(abs(A),2);                               

                % Signal Smoothing
                coef = ones(1,avg_value_set(iii))./avg_value_set(iii);                      % FIR filter for signal averaging
                avg_estimate = filter(coef,1,estimate);                                     % apply filter to signal

                % FIR filter delay
                delay = round((length(coef)-1)/2,0);                                        % Estimate introduced delay
                estimate = estimate(1:end-delay);                                           % Take into account FIR filter delay
                avg_estimate(1:delay) = [];
                
                % Estimation extraction
                [~,peak_index] = max(estimate);
                [~,peak_index_avg] = max(avg_estimate);

                % Error calculation
                if (peak_index > (T_symb+CP)/2)
                    peak_index = peak_index - (T_symb+CP);
                end

                if (peak_index_avg > (T_symb+CP)/2)
                    peak_index_avg = peak_index_avg - (T_symb+CP);
                end
                
                E(i,ii) = n_tram - peak_index;
                E((i+iii*length(snr_set)),ii) = n_tram - peak_index_avg;
                
            end
        end
            
%             Plot Worst and Best case (only used to check signal averaging delay)
%             if i==1 && ii==1
%                 figure;
%                 hold on;
%                 plot(estimate);
%                 plot(avg_estimate);
%                 title(['Worst case estimate. SNR = ',num2str(snr_set(i))],['Delay is = ',num2str(n_tram)]);
%                 legend('Original',num2str(avg_value));
%             end
%             
%             if i==length(snr_set) && ii==1
%                 figure;
%                 hold on;
%                 plot(estimate);
%                 plot(avg_estimate);
%                 title(['Best case estimate. SNR = ',num2str(snr_set(i))],['Delay is = ',num2str(n_tram)]);
%                 legend('Original',num2str(avg_value));
%             end
end

% RMS calculation (column dimension due to how data has been allocated)
for i=1:(length(avg_value_set)+1)
    RMS_E(((i-1)*length(snr_set) + 1):i*length(snr_set),:) = rms(E(((i-1)*length(snr_set) + 1):i*length(snr_set),:),2);
end

% RMS error plot
figure;
hold on;
grid on;

for ii=1:(length(avg_value_set)+1)
    plot(snr_set,RMS_E(((ii-1)*length(snr_set) + 1):ii*length(snr_set),:));
end
    
yline(max_accepted_rmse,'--','Threshold','LineWidth',2);
legend(['No Averaging';append('Averaging ',cellstr(int2str(avg_value_set')))]);
title(['RMS Error (rational delay). Number of simulations per SNR value = ',num2str(n_simul)]);
ylabel('Sample Error');
xlabel('SNR');

toc
