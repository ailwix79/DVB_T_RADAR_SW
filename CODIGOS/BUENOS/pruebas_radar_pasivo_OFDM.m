clear all;
close all;
clc;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% OFDM signal generation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

T_symb  = 8192;
CP      = T_symb/32;
n_symb  = 128; % Number of symbols

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% M-QAM symbols over each carrier
M       = 64;
symbols = (randi(sqrt(M),[T_symb,n_symb])-(sqrt(M)+1)/2) + 1j*(randi(sqrt(M),[T_symb,n_symb])-(sqrt(M)+1)/2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Add continual pilot information on pseudorandom carrier locations indexes
pilot_data = randi(2,[256,1])-1.5;
pilot_carr_aux = unique(randi(T_symb,[1,512]),'stable');
pilot_carr = pilot_carr_aux(logical((pilot_carr_aux<(T_symb/2-round(T_symb*0.05)))+(pilot_carr_aux>(T_symb/2+round(T_symb*0.05)))));
pilot_carr = pilot_carr(1:256);
symbols(pilot_carr,:) = repmat(pilot_data,1,n_symb);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% High frequency symbols are deteled
symbols(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Perform IFFT transform and add ciclic prefix
s_tx = ifft(symbols,T_symb,1);
s_tx = [s_tx(end-CP+1:end,:); s_tx];
s_tx = s_tx(:);
s_tx = s_tx/sqrt(mean(abs(s_tx).^2)); % Unit power normalization

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Channel Model
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Direct Path to noise power ratio
SNR = 40;
snr = 10^(SNR/10);

% Set antenna transmitter position
r_tx = [-5e3;0;0];

% Set passive radar receiver position
r_rx = [5e3;0;0];

% Cinematic profile of the target and initial position
vo = [0;-300;0]; % m/s
ro = [0;10e3;3e3]; % m

P = 1; % Simulation step (samples)
fs = 10e6; % Hz
f = 300e6;

s_tx = [s_tx; zeros(P-mod(length(s_tx),P),1)];

s_rx = zeros(length(s_tx),1);
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
s_rx = s_tx + 0.001*s_rx;
% h = 0.1.^(0:4095);
% s_rx = filter(h(:),1,s_tx) + 0.0001*s_rx;

s_rx = s_rx + (1/sqrt(snr))*(randn(size(s_rx)) + 1j*randn(size(s_rx)))/sqrt(2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Signal Separation Algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cross-Correlation of the received signal
s_x = conj([s_rx; zeros(T_symb,1)]).*[zeros(T_symb,1); s_rx];
s_x = conv(s_x, ones(CP,1));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sliding window peak detector (Coarse Symbol Synchronization)
W = T_symb+CP; % Sliding window length
sof = []; % Start of Frame index vector
amp = []; % Start of Frame Correlation amplitude

x_buff = s_x(1:min(T_symb+CP,length(s_x)));
[~,r_idx] = max(abs(x_buff));
val = x_buff(r_idx);
sof(end+1) = r_idx;
amp(end+1) = val;

while sof(end)+(T_symb+CP)/2 < length(s_x)
    x_buff = s_x(sof(end)+(T_symb+CP)/2:min((sof(end)+(T_symb+CP)/2 + T_symb+CP - 1),length(s_x)));
    [~,r_idx] = max(abs(x_buff));
    val = x_buff(r_idx);
    sof(end+1) = r_idx+sof(end)+(T_symb+CP)/2-1;
    amp(end+1) = val;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Perfect Timing synchronization 
% sof = ((1:n_symb)-1)*(T_symb+CP);
% amp = s_x(sof+1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Data demodulation and reconstruction

% Intermediate variables definition
y           = complex(zeros(T_symb,1));
Y           = complex(zeros(T_symb,1));
H           = complex(zeros(T_symb,1));
Yeq         = complex(zeros(T_symb,1));
Ydec        = complex(zeros(T_symb,1));
Yobs        = complex(zeros(T_symb,1));
y_rx        = complex(zeros(T_symb,1));
y_rx_cp     = complex(zeros(T_symb+CP,1));
y_dp        = complex(zeros(length(s_rx)+(T_symb+CP),1)); % reconstructed Reference channel signal
y_obs       = complex(zeros(T_symb,1));
y_obs_cp    = complex(zeros(T_symb+CP,1));
y_oc        = complex(zeros(length(s_rx)+(T_symb+CP),1)); % reconstructed Observation channel signal

for i = 1:(length(sof)-1)
    H(:) = complex(zeros(T_symb,1));
    % 1) Remove cyclic prefix and Perform FFT
    y(:) = zeros(T_symb,1);
    y(1:min(T_symb,length(s_rx)-(sof(i)+CP)),:) = s_rx(sof(i)+CP+1:min(sof(i)+CP+T_symb,length(s_rx)),:);
    Y(:) = fft(y,T_symb);
    % 2) Calculate equalizer coefficients
    H(pilot_carr,:) = Y(pilot_carr,:)./pilot_data;
%     H(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 1;
    H(:) = interp1(sort(pilot_carr),H(sort(pilot_carr),:),1:T_symb,'linear','extrap');
    % 3) Fine Symbol Synchronization (from the impulse response of H)
    [~, tau] = max(abs(ifft(H,T_symb)));
    if tau >= T_symb/2
        tau = tau-T_symb;
    end
    sof(i) = round(sof(i)+tau-1);
%     tau = mean(-1*diff(unwrap(angle(H(1:round(T_symb*0.25))))*T_symb/(2*pi))); 
%     sof(i) = round(sof(i)+tau);
    y(:) = zeros(T_symb,1);
    y(1:min(T_symb,length(s_rx)-(sof(i)+CP)),:) = s_rx(sof(i)+CP+1:min(sof(i)+CP+T_symb,length(s_rx)),:);
    Y(:) = fft(y,T_symb);
    % 5) Re-Calculate equalizer coefficients
    H(:) = complex(zeros(T_symb,1));
    H(pilot_carr,:) = Y(pilot_carr,:)./pilot_data;
%     H(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 1;
    H(:) = interp1(sort(pilot_carr),H(sort(pilot_carr),:),1:T_symb,'linear','extrap');
    % 6) Zero-forcing equalizer
    Haux = H;
    Haux(H==0) = 1;
    Yeq(:) = Y./Haux;
    % 7) M-QAM Hard Decision Decoder
    Ydec(:) = round(Yeq + (sqrt(M)/2-0.5)*(1+1j));
    Ydec(real(Ydec)>(sqrt(M)-1)) = (sqrt(M)-1) + 1j*imag(Ydec(real(Ydec)>(sqrt(M)-1)));
    Ydec(real(Ydec)<0) = 1j*imag(Ydec(real(Ydec)<0));
    Ydec(imag(Ydec)>(sqrt(M)-1)) = 1j*(sqrt(M)-1) + real(Ydec(imag(Ydec)>(sqrt(M)-1)));
    Ydec(imag(Ydec)<0) = real(Ydec(imag(Ydec)<0));
    Ydec(:) = Ydec - (sqrt(M)/2-0.5)*(1+1j);
    Ydec(pilot_carr) = pilot_data; 
    Ydec(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 0;
    % 8) Signal reconstruction
    y_rx(:) = ifft(Ydec);
    % 9) Add cyclic prefix
    y_rx_cp(:) = [y_rx(end-CP+1:end); y_rx];
    % 10) Save reconstructed signal
    y_dp(sof(i)+1:sof(i)+CP+T_symb) = y_rx_cp;
    % 11) Obtain Clutter-free observation channel signal
    Yobs(:) = Y-Ydec.*Haux;
    y_obs(:) = ifft(Yobs);
    y_obs_cp(:) = [y_obs(end-CP+1:end); y_obs];
    y_oc(sof(i)+1:sof(i)+CP+T_symb) = y_obs_cp;
end

y_oc(isnan(y_oc))=0;
y_dp(isnan(y_dp))=0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Correlation processing (Batch Algorithm)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

L = 2^12; % Decimation factor (Batch block processing length)

s_rx = y_oc;
s_tx = y_dp(1:min(length(y_dp),length(s_rx)));
s_tx = y_dp;

s_rx = [s_rx(:); zeros(L-mod(length(s_rx),L),1)];
s_tx = [s_tx(:); zeros(length(s_rx)-length(s_tx),1)];

s_rx = reshape(s_rx,L,length(s_rx)/L);
s_tx = reshape(s_tx,L,length(s_tx)/L);

cols = size(s_rx,2); 

Y = zeros(L,cols);

Lfft_aux = 2^(nextpow2(2*L-1));
y_aux = zeros(Lfft_aux,1);

% Time-domain processing
for i = 1:cols
    y_aux(:) = ifft(fft(s_rx(:,i),Lfft_aux).*conj(fft(s_tx(:,i),Lfft_aux)),Lfft_aux);
    % y_aux(:) = xcorr(s_rx(:,i),s_tx(:,i));
    % Y(:,i) = y_aux(L:2*L-1);
    Y(:,i) = y_aux(1:L);
end

% Frequency Domain processing
K = 2^(nextpow2(cols)); % FFT length
Y = [Y, zeros(L,K-cols)]; % zero-padding

% Window function applied to reduce the sidelobe level
w = chebwin(cols,80);
w = w(:)/sum(w);
w = [w; zeros(K-cols,1)];

% Doppler processing
Y = Y.*repmat(w.',L,1);
Y = fftshift(fft(Y,K,2),2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% PLOTs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f1 = figure();
a1 = axes();
surf(a1,((0:(K-1))/K-0.5)*(fs/L)*3e8/(2*f),(0:(L-1))*3e8/(2*fs),10*log10((abs(Y)).^2))
xlim(a1,[-3e8*fs/(4*L*f),3e8*fs/(4*L*f)])
ylim(a1,[0,(L-1)*3e8/(2*fs)])
ylabel(a1,'Bistatic Range (m)','Interpreter','Latex')
xlabel(a1,'Bistatic Velocity (m/s)','Interpreter','Latex')
zlabel(a1,'Power Level (dB)','Interpreter','Latex')
title(a1,['DNR = ',num2str(SNR),' dB, $T_{coh} = ',num2str(cols*L/fs,'%4.2f'),'$ sec.'],'Interpreter','Latex')
shading interp
colormap jet