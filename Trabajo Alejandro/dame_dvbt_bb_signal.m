% DAME_DVBT_BB_SIGNAL
%
% y = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T)
%
% @Brief: This function generates a dvb-t signal vector according to the
% DVB-T standard definition (ETSI EN 300 744 V1.6.2 (2015-10)).
%
% Input parameters:
%   - BW            : Channel bandwidth in MHz (8, 7, 6 or 5)
%   - tx_mode       : '2K', '4K' or '8K' transmitter mode
%   - frame_offset  : frame offset within a superframe structure
%   - guard         : Guard period (as a fraction of the symbol period)
%   - mod_type      : Data modulation format ('QPSK', '16-QAM' or '64-QAM')
%   - alpha         : Alpha constant (see ETSI EN 300 744 V1.6.2) (1, 2 or 4)
%   - T             : Vector length (in seconds)
%
% Output parameters:
%   - y             : output dvb-t signal vector
% 
% @Author: Carlos García de la Cueva
% @Revision:    v1.0    -   24/05/2022
%               Progress update rules:
%               v1.xx   -   Minor changes related to release #1
%               v2.00   -   Mayor release update
%

% ALEJANDRO
% Changed function output to [y,T_symb] for cyclic guard length calculation in
% main code.

% TEST VALUES
% BW = 8;                                 % Channel selection 5,6,7,8 (Mhz)
% T = 1e-5;                               % Symbol duration (1/fs)
% tx_mode = '8K';                         % Transmitter mode '2K','4K','8K' 
%                                         % (described in the standart, different number of carriers and Useful Time values for each mode)
% frame_offset = T/16;                    % Frame Offset within a superframe (PREGUNTAR)
% guard = 1/32;                           % Guard interval length (fraction of T, page 33 table 14)
% mod_type = '64-QAM';                    % Symbol modulation 'QPSK','16-QAM','64-QAM'
% alpha = 2;                              % Normalization factor for the modulation of the OFDM symbol (page 27)

function [y,Tu] = dame_dvbt_bb_signal(BW, tx_mode, frame_offset, guard, mod_type, alpha, T)

% DVB-T Transmitter Sampling Frequency
fs = [64/7, 8, 48/7, 40/7, 5]*1e6; % Hz for an 8, 7, 6 and 5 MHz channel respectively
switch BW
	case 8
		fs = fs(1);
	case 7
		fs = fs(2);
	case 6
		fs = fs(3);
	case 5
		fs = fs(4);
end

% DVB-T Modulation type and normalization factor
% Pagina 27 normalization factor for data symbols

switch mod_type
	case 'QPSK'
		M = 4;
		c = 1/sqrt(2);
	case '16-QAM'
		M = 16;
		switch alpha
			case 1
				c = 1/sqrt(10);
			case 2
				c = 1/sqrt(20);
			case 4
				c = 1/sqrt(52);
		end
	case '64-QAM'
		M = 64;
		switch alpha
			case 1
				c = 1/sqrt(42);
			case 2
				c = 1/sqrt(60);
			case 4
				c = 1/sqrt(108);
		end
end

% Overall signal simulation period (expressed in samples)
T = T*fs;

% Continual pilot carrier indexes
% Pagina 29 carrier indices for continual pilot carriers

cont_pilot_carriers = [ 0       48      54      87      141     156     192 ...
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

% transmission parameter signalling carriers
% Pagina 30 TPS carrier indices
% TPS: 68 OFDM symbols are one OFDM frame, four frames are a super frame
% Each OFDM symbols conveys one TPS bit, each TPS block has 68 bits
%   1 initialization bit
%   16 synchronization bits
%   37 information bits
%   14 redundancy bits
%   remaining 6 bits set to zero

% 
tps_carriers 		= [	34      50      209     346     413     569     595  ...
						688 	790     901     1073    1219    1262    1286 ...
						1469    1594 	1687    1738    1754    1913    2050 ...
						2117    2273    2299 	2392    2494    2605    2777 ...
						2923    2966    2990    3173 	3298    3391    3442 ...
						3458    3617    3754    3821    3977 	4003    4096 ...
						4198    4309    4481    4627    4670    4694 	4877 ...
						5002    5095    5146    5162    5321    5458    5525 ...
						5681    5707    5800    5902    6013    6185    6331 ...   
						6374  	6398    6581    6706    6799];

% Kmax: maximum carrier value
% Tu: useful part duration, page 27
% Ts: total transmission time
% Tf: individual frame duration
% The spacing between carriers is 1/Tu while the spacing between Kmin and
% Kmax is (K-1)/Tu
% Guard interval: duration of guard interval / Tu

% PREGUNTAR L_inf
% L_inf numero de portadoras que contienen las portadoras que contienen
% informacion
switch tx_mode
    case '2K'
        % 2K Mode
        cont_pilot_carriers = cont_pilot_carriers(1:45);
        tps_carriers = tps_carriers(1:17);
        Kmax = 1704;    % Tambien valdría max(cont_pilot_carriers)
		Tu = 2048;
		L_inf = 1512;
    case '4K'
        % 4K Mode
        cont_pilot_carriers = cont_pilot_carriers(1:89);
        tps_carriers = tps_carriers(1:34);
        Kmax = 3408;    % Tambien valdría max(cont_pilot_carriers)
		Tu = 4096;
    case '8K'
        % 8K Mode
        Kmax = 6816;
		Tu = 8192;
		L_inf = 6048;
    otherwise
        error('Invalid MODE identifier!') 
end

% Available data carriers indexes
% Generate negative and positive carriers
% Example Kmax=6816 data_carrier_indexes = [-3408 , 3408]
data_carriers_indexes = (0:Kmax)-Kmax/2;

% PREGUNTAR
% evitar fftshift
data_carriers_indexes(data_carriers_indexes<0) = data_carriers_indexes(data_carriers_indexes<0) + Tu;

% Number of symbols to be generated

% PREGUNTAR
% calculo del número de símbolos (Tiempo total entre tiempo de transmision
% útil
nsymb = ceil(T/(Tu*(1+guard)));

% Initialize information matrix
Y = complex(zeros(Tu,nsymb));

% Pseudorandom binary sequence
% x^11 + x^2 + 1
% Generation schema on page 28
wk = prbs_dvbt(Kmax+1);

% Get information of continual pilots
% Continual and scattered pilots are modulated according to a PRBS sequence
% It also governs the starting phase of the TPS information.
cont_pilot_inf = wk(cont_pilot_carriers+1);

% Transform pilot carrier indexes to FFT indexes
cont_pilot_carriers(:) = cont_pilot_carriers - Kmax/2;
cont_pilot_carriers(cont_pilot_carriers<0) = cont_pilot_carriers(cont_pilot_carriers<0) + Tu;

% Add continual pilots to the information matrix at boosted power level
% Formula en pagina 28, situate scattered pilot cells
Y(cont_pilot_carriers+1,:) = (4/3)*2*(0.5-repmat(cont_pilot_inf(:),1,nsymb));

tps_carriers(:) = tps_carriers - Kmax/2;
tps_carriers(tps_carriers<0) = tps_carriers(tps_carriers<0) + Tu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter Pilots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 0:floor(Kmax/12);
tmp = zeros(1,length(p));
scatt_pilot_carriers = [12*p(1:end-1);3+12*p(1:end-1);6+12*p(1:end-1);9+12*p(1:end-1)];
scatt_pilot_carriers_aux = scatt_pilot_carriers - Kmax/2;
scatt_pilot_carriers_aux(scatt_pilot_carriers_aux<0) = scatt_pilot_carriers_aux(scatt_pilot_carriers_aux<0) + Tu;

scatt_pilot_inf = zeros(1,size(scatt_pilot_carriers,2));

data_carriers = zeros(Tu,1);
data = complex(zeros(L_inf,1));

% Scatter pilots for every frame within a super frame
for l = frame_offset:(nsymb+frame_offset-1)
    tmp(:) = 3*mod(l,4)+12*p;
	scatt_pilot_inf(:) = wk(scatt_pilot_carriers(mod(l,4)+1,:)+1);
	Y(scatt_pilot_carriers_aux(mod(l,4)+1,:)+1,l-frame_offset+1) = (4/3)*2*(0.5-scatt_pilot_inf(:));
	% Get data available carriers
	data_carriers(:) = zeros(Tu,1);
    data_carriers(data_carriers_indexes+1) = 1;
	data_carriers(cont_pilot_carriers+1) = 0;
	data_carriers(scatt_pilot_carriers_aux(mod(l,4)+1,:)+1) = 0;
	data_carriers(tps_carriers+1) = 0;
	data(:) = 2*((randi(sqrt(M),[L_inf,1])-(sqrt(M)+1)/2) + 1j*(randi(sqrt(M),[L_inf,1])-(sqrt(M)+1)/2));
	data(:) = data + (alpha-1)*(sign(real(data))+1j*sign(imag(data)));
	Y(logical(data_carriers),l-frame_offset+1) = data*c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cyclic Prefix insertion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = ifft(Y);
y = [y(end-Tu*guard+1:end,:);y];
y = y(:);