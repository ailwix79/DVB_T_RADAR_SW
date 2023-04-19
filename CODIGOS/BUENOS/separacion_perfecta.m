
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % DVB_T Signal Generation
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    snr_db = 60;
    BW = 8;                                 % Channel selection 5,6,7,8 (Mhz)
    n_symb = 128;                           % Number of symbols
    fs = (64/7)*1e6;                        % Sampling frequency (Hz)
    tx_mode = '8K';                         % Transmitter mode '2K','4K','8K', different number of carriers and Tu values for each mode...
    frame_offset = 0;                       % Frame offset
    guard = 1/32;                           % Guard interval length (fraction of T, page 33 table 14)
    T_symb = 8192;                          % Useful symbol time in samples
    CP = guard*T_symb;                      % Guard length in samples
    mod_type = '64-QAM';                    % Symbol modulation 'QPSK','16-QAM','64-QAM'
    M = 64;
    alpha = 2;                              % Normalization factor for the modulation of the OFDM symbol (page 27)
    T = (n_symb*(T_symb + CP))/fs;          % Signal duration (in seconds)

    [s_tx,~,pilot_carr,data] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Channel Model
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Direct Path to noise power ratio
    snr = 10^(snr_db/10);

    % Set antenna transmitter position
    r_tx = [-5e3;0;0];

    % Set passive radar receiver position
    r_rx = [5e3;0;0];

    % Cinematic profile of the target and initial position
    vo = [200;0;0]; % m/s
    ro = [5e3;0;3e3]; % m

    P = 1; % Simulation step (samples)
    fs = (64/7)*1e6; % Hz
    f = 300e6;

    s_tx = [s_tx; zeros(P-mod(length(s_tx),P),1)];

    s_rx = zeros(length(s_tx),1); % generate space for the incoming signal
    s_aux = zeros(P,1);

    r = 0;
    dt = 0;

    for i = 1:(length(s_tx)/P)
        r(:) = sqrt(sum((r_tx-(ro + vo*(i-1)*P/fs)).^2))+sqrt(sum((r_rx-(ro + vo*(i-1)*P/fs)).^2));
        dt(:) = fs*2*r/3e8;
        if (ceil(dt)+((i-1)*P+1)) <= length(s_rx)
            s_aux(:) = s_tx(((i-1)*P+1):i*P);
            s_rx(ceil(dt)+((i-1)*P+1):min(ceil(dt)+i*P,length(s_rx))) = s_rx(ceil(dt)+((i-1)*P+1):min(ceil(dt)+i*P,length(s_rx))) + s_aux(1:min(ceil(dt)+i*P,length(s_rx))-(ceil(dt)+((i-1)*P+1))+1)*exp(-1j*2*pi*f*r/3e8);
        end
        % Print progress update
        if floor(100*i/(length(s_tx)/P)) > floor(100*(i-1)/(length(s_tx)/P))
          clc;
          fprintf('Generando el escenario ... (%3.0f %%)\n',floor(100*i/(length(s_tx)/P)));
        end
    end

    % s_rx = s_tx;
    % h = 0.1.^(0:4095);
    % s_rx = filter(h(:),1,s_tx) + 0.0001*s_rx;

    s_rx = s_rx + (1/sqrt(snr))*(randn(size(s_rx)) + 1j*randn(size(s_rx)))/sqrt(2);
    
    
% BATCH ALGORITHM
L = 2^12;
lfft = 2^(nextpow2(2*L-1));
factor_interp = 1;
[Y2,cols] = batch_corr_algorithm(s_rx,s_tx,L,factor_interp);

% RESULTS
figure('Name','Alejandro');
a1 = axes();
surf(a1,((0:(2^nextpow2(cols)-1))/2^nextpow2(cols)-0.5)*(fs/L)*3e8/(2*f),(0:(L-1))*3e8/(2*fs),10*log10((abs(Y2)).^2))
xlim(a1,[-3e8*fs/(4*L*f),3e8*fs/(4*L*f)])
ylim(a1,[0,(L-1)*3e8/(2*fs)])
ylabel(a1,'Bistatic Range (m)','Interpreter','Latex')
xlabel(a1,'Bistatic Velocity (m/s)','Interpreter','Latex')
zlabel(a1,'Power Level (dB)','Interpreter','Latex')
title(a1,['SNR = ',num2str(snr_db),' dB, $T_{coh} = ',num2str(cols*L/fs,'%4.2f'),'$ sec.'],'Interpreter','Latex')
shading interp
colormap jet