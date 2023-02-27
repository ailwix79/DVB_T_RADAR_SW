% Channel estimation and signal separation
% Scattered and continual pilots declaration
% DVB-T signal mode is 8K

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

[s_rx,~] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T);

% Continual pilots (these are fixed by the standart)
% In 8K all pilots are used
pilot =               [ 0       48      54      87      141     156     192 ...
                        201     255     279     282     333     432     450 ...
                        483     525     531     618     636     714     759 ...
                        765     780     804     873     888     918     939 ...
                        942     969     984     1050    1101    1107    1110 ...
                        1137    1140    1146    1206    1269    1323    1377 ...
                        1491    1683    1704    1752    1758    1791    1845 ... 
                        1860    1896    1905    1959    1983    1986    2037 ...
                        2136    2154    2187    2229    2235    2322    2340 ...
                        2418    2463    2469    2484    2508    2577    2592 ...
                        2622    2643    2646    2673    2688    2754    2805 ...
                        2811    2814    2841    2844    2850    2910    2973 ...
                        3027    3081    3195    3387    3408    3456    3462 ...
                        3495    3549    3564    3600    3609    3663    3687 ...
                        3690    3741    3840    3858    3891    3933    3939 ...
                        4026    4044    4122    4167    4173    4188    4212 ...
                        4281    4296    4326    4347    4350    4377    4392 ...
                        4458    4509    4515    4518    4545    4548    4554 ...
                        4614    4677    4731    4785    4899    5091    5112 ...
                        5160    5166    5199    5253    5268    5304    5313 ...
                        5367    5391    5394    5445    5544    5562    5595 ...
                        5637    5643    5730    5748    5826    5871    5877 ...
                        5892    5916    5985    6000    6030    6051    6054 ...
                        6081    6096    6162    6213    6219    6222    6249 ...
                        6252    6258    6318    6381    6435    6489    6603 ...
                        6795    6816];
                    
% Scatter pilots
% - p: scatter pilot index vector, must be all possible values greater than or equal to zero, provided that the resulting
%      value for k does not exceed the valid range [K min ; K max].
% - k: scatter pilot index subset (frequency bins). k = Kmin + 3 × (l mod 4) + 12p
% - l: index symbol, ranging from 0 to 67
% - v: expected "boosted" modulation value of scatter pilot

% Maximum k value will be kmax = 0 + 3 * mod(67,4) + 12p => kmax = 9 + 12p
% possible mod(l,4) values are 0,1,2,3 (each superframe)
% therefore all possible k formulas are

% - kmin = 12p          Sub frame 1 mod(l,4) = 0
% - k = 3 + 12p         Sub frame 2 mod(l,4) = 1
% - k = 6 + 12p         Sub frame 3 mod(l,4) = 2
% - kmax = 9 + 12p      Sub frame 4 mod(l,4) = 3

% Hence, p must be a vector that makes sure k values do not leave the range
% [Kmin, Kmax] that is to say that kmin >= Kmin and kmax <= Kmax
% 12p >= Kmin and 9 + 12p <= Kmax
% therefore p_min = 0 and p_max = (Kmax-9)/12

Kmin = 0;
Kmax = 6816;

l = 0:67;
p = 0:ceil((Kmax-9)/12);
k = [(12*p)' , (3 + 12*p)' , (6 + 12*p)' , (9 + 12*p)'];

wk = zeros(1,Kmax+1);

mem = (ones(1,11));

for i = 1:(Kmax+1)
    wk(i) = mem(1);
    mem(:) = [mem(2:end) bitxor(mem(1),mem(3))];
end

v = unique((4/3)*2*(0.5-wk));

s_x = conj([s_rx; zeros(T_symb,1)]).*[zeros(T_symb,1); s_rx];      
s_x = conv(s_x, ones(CP,1)); 

[pks,locs] = findpeaks(s_x,'MinPeakDistance',(T_symb+CP)/2);

locs(pks < (mean(pks)/2)) = [];
pks(pks < (mean(pks)/2)) = [];

for i=1:((length(s_rx)/(T_symb + CP))-1)
    
    if (mod(i-1,N) == 0)
        
        H_full = zeros(T_symb,1);

        for n=1:N
            % Eliminate CP
            num_sample = i-1+n;
            symb = zeros(T_symb,1);
            symb(1:min(T_symb,length(s_rx)-(locs(num_sample)+CP)),:) = s_rx((locs(num_sample)+CP+1):min(locs(num_sample)+CP+T_symb,length(s_rx)),:);

            % Obtain amplitude of subcarriers
            C = symb(pilot_carr);

            % Apply FFT to OFDM symbol
            X = fft(symb, T_symb);
            
            % Multiplicacion en puntos de pilotos dispersos, donde espero
            % recibir el simbolo BPSK
            
            % TRABAJAR AQUI
            % 4 VECTORES PARA CADA UNO DE LOS POSIBLES ÍNDICES DISPERSOS
            % Correlacion
            % Deteccion de Patrones
            
            % Estimate channel responde (FFT of sub carrier divided by amplitude of subcarrier)
            H = X(pilot_carr + indices_dispersas)./C;

            % Interpolate with pilot subcarrier indexes to obtain full channel response
            H_full = H_full + interp1(pilot_carr,H,1:T_symb,'linear','extrap')';

        end
        
        % Obtain averaged channel estimate for N symbols
        H_average = H_full./N;
    end
    
    % Equalization
    Yeq = X./H_average;
    
    % TODO
    % - Decisor QAM
    % - Xf = X - H*Cq
    % - Inserción CP
    % - IFFT
    
end
