clear all;
clc;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% OFDM signal generation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

T_symb  = 8192;                 % Duración de un símbolo, longitud de la IFFT, longitud de un simbolo OFDM sin tener en cuenta CP
CP      = T_symb/32;            % longitud prefijo cíclico, longitud sym OFDM total T_symb + CP
n_symb  = 128;                  % numero de símbolos OFDM o portadoras piloto

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% M-QAM symbols over each carrier
M       = 16; % Numero de puntos en la constelación
symbols = (randi(sqrt(M),[T_symb,n_symb])-(sqrt(M)+1)/2) + 1j*(randi(sqrt(M),[T_symb,n_symb])-(sqrt(M)+1)/2); % Constelación QAM

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Add continual pilot information on pseudorandom carrier locations indexes


pilot_data = randi(2,[256,1])-1.5;                                      % Datos a enviar
pilot_carr_aux = unique(randi(T_symb,[1,512]),'stable');                % Portadoras (piloto y sub)

% PREGUNTAR
pilot_carr = pilot_carr_aux(logical((pilot_carr_aux<(T_symb/2-round(T_symb*0.05)))+(pilot_carr_aux>(T_symb/2+round(T_symb*0.05))))); % pilot carriers

pilot_carr = pilot_carr(1:256); 	% sub carriers
symbols(pilot_carr,:) = repmat(pilot_data,1,n_symb); %data modulation M-QAM, each carrier (index) for each data element. The resulting constellation is the transmission

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% High frequency symbols are deleted
symbols(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 0;     % We are only interested in the lowband signal

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Perform IFFT transform and add ciclic prefix
s_tx = ifft(symbols,T_symb,1);
s_tx = [s_tx(end-CP+1:end,:); s_tx];    % cyclic prefix used for protectiong against intersymbol error
s_tx = s_tx(:);
s_tx = s_tx/sqrt(mean(abs(s_tx).^2)); % Unit power normalization

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure;
plot(10*log10(abs(s_tx)));
title('Señal generada OFDM');
% Channel Model
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Direct Path to noise power ratio
SNR = 60;
snr = 10^(SNR/10);

% Set antenna transmitter position
r_tx = [-5e3;0;0];

% Set passive radar receiver position
r_rx = [5e3;0;0];

% Cinematic profile of the target and initial position
vo = [0;200;0]; % m/s
ro = [0;10e3;3e3]; % m

P = 1; % Simulation step (samples)
fs = 10e6; % Hz
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
s_rx = s_tx + 0.001*s_rx;
% h = 0.1.^(0:4095);
% s_rx = filter(h(:),1,s_tx) + 0.0001*s_rx;

s_rx = s_rx + (1/sqrt(snr))*(randn(size(s_rx)) + 1j*randn(size(s_rx)))/sqrt(2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Signal Separation Algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % Cambiar sincronismo fino, sincronismo grueso, estimación de canal
% 
% % Cross-Correlation of the received signal
% s_x = conj([s_rx; zeros(T_symb,1)]).*[zeros(T_symb,1); s_rx];      
% s_x = conv(s_x, ones(CP,1)); 
% 
% figure;
% plot(abs(s_x));
% title('Cross correlation');
% xlabel('Samples');
% 
% % Detect cross correlation peaks (pilot carrier position)
% % Perfect Timing synchronization 
% [~,ind] = max(1:(T_symb + CP));                      % find first pilot carrier index
% I = ind:(T_symb + CP):ind*n_symb; I = I';            % all pilots are separated T_symb + CP
% V = s_rx(I);                                         % get pilot amplitude
% 
% % CFO correction should occur here, before the actual CP removal and
% % demodulation
% 
% for i=1:((length(s_x)/(T_symb + CP))-1)
%     % Eliminate cyclic prefix
%     y(1:min(T_symb,length(s_rx)-(I(i)+CP)),:) = s_rx(I(i)+CP+1:min(I(i)+CP+T_symb,length(s_rx)),:);
%     
% end

s_x = conj([s_rx; zeros(T_symb,1)]).*[zeros(T_symb,1); s_rx];      
s_x = conv(s_x, ones(CP,1)); 

plot(abs(s_x));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sliding window peak detector (Coarse Symbol Synchronization)
% Ventana deslizante que detecta el máximo de la señal de la parte que
% filtra
% Ventana grande: picos sin detectar
% Ventana pequeña: pseudo picos

W = T_symb+CP; % Longitud de la ventana, misma longitud que un frame de OFDM (longitud simbolo + CP) 129 ventanas (length(s_x)/(T_symb + CP))
sof = []; % Comienzo ventana
amp = []; % Final ventana

x_buff = s_x(1:min(T_symb+CP,length(s_x)));         % buffer
[~,r_idx] = max(abs(x_buff));                       % indice del mayor valor del buffer
val = x_buff(r_idx);                                % valor máximo del buffer.
sof(end+1) = r_idx;                                 % inicializar inicio ventana
amp(end+1) = val;                                   % inicializar final ventana

% Este bucle while va moviendo la ventana por la señal

while sof(end)+(T_symb+CP)/2 < length(s_x)
    % El buffer va tomando muestras de la señal
    x_buff = s_x(sof(end)+(T_symb+CP)/2:min((sof(end)+(T_symb+CP)/2 + T_symb+CP - 1),length(s_x)));
    [~,r_idx] = max(abs(x_buff));
    val = x_buff(r_idx);
    sof(end+1) = r_idx+sof(end)+(T_symb+CP)/2-1;
    amp(end+1) = val;
    
    % El final de la ventana se pone en el maximo valor detectado, el
    % principio en su índice
end

% sof me dice donde están los valores máximos y amp me dice que valor
% tienen. El proposito de la parte anterior es detectar las portadoras
% piloto (tonos empleados para la equalización más adelante)

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

% para cada símbolo OFDM de longitud T_symb

for i = 1:(length(sof)-1)   %también valdria i = 1:(length(s_x)/(T_symb + CP))
    H(:) = complex(zeros(T_symb,1));
    % 1) Remove cyclic prefix and Perform FFT
    y(:) = zeros(T_symb,1);
    y(1:min(T_symb,length(s_rx)-(sof(i)+CP)),:) = s_rx(sof(i)+CP+1:min(sof(i)+CP+T_symb,length(s_rx)),:);      % elimina el CP
    Y(:) = fft(y,T_symb);
    % 2) Calculate equalizer coefficients El canal puede introducir desfase
    % en las subcarriers. Las portadoras piloto me dicen cuanto tengo que
    % corregir.
    
    H(pilot_carr,:) = Y(pilot_carr,:)./pilot_data;             % se estiman los valores de equalización con los valores de las portadoras piloto
%     H(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 1;
    H(:) = interp1(sort(pilot_carr),H(sort(pilot_carr),:),1:T_symb,'linear','extrap');          % H es la respuesta al impulso del canal
    % 3) Fine Symbol Synchronization (from the impulse response of H)
    
    % PREGUNTAR correción de fase
    [~, tau] = max(abs(ifft(H,T_symb)));
    if tau >= T_symb/2
        tau = tau-T_symb;
    end
    sof(i) = round(sof(i)+tau-1);
%     tau = mean(-1*diff(unwrap(angle(H(1:round(T_symb*0.25))))*T_symb/(2*pi))); 
%     sof(i) = round(sof(i)+tau);

    % se ha modificado sof (posiciones de las portadoras piloto) por la
    % sincronización. Es necesario recalcular donde estran los CP para
    % eliminarlos
    
    
    %=======
    
    
    y(:) = zeros(T_symb,1);
    y(1:min(T_symb,length(s_rx)-(sof(i)+CP)),:) = s_rx(sof(i)+CP+1:min(sof(i)+CP+T_symb,length(s_rx)),:);
    Y(:) = fft(y,T_symb);
    
    % 5) Re-Calculate equalizer coefficients Con la sincronización
    % realizada se recalculan los parámetros de equalización
    
    H(:) = complex(zeros(T_symb,1));
    H(pilot_carr,:) = Y(pilot_carr,:)./pilot_data;          % Channel estimation LS method
%     H(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 1;
    H(:) = interp1(sort(pilot_carr),H(sort(pilot_carr),:),1:T_symb,'linear','extrap');
    
    % 6) Zero-forcing equalizer
    Haux = H;
    Haux(H==0) = 1;         %PREGUNTAR
    Yeq(:) = Y./Haux;
    
    % 7) M-QAM Hard Decision Decoder Decodificacion por decisor QAM
    Ydec(:) = round(Yeq + (sqrt(M)/2-0.5)*(1+1j));
    Ydec(real(Ydec)>(sqrt(M)-1)) = (sqrt(M)-1) + 1j*imag(Ydec(real(Ydec)>(sqrt(M)-1)));
    Ydec(real(Ydec)<0) = 1j*imag(Ydec(real(Ydec)<0));
    Ydec(imag(Ydec)>(sqrt(M)-1)) = 1j*(sqrt(M)-1) + real(Ydec(imag(Ydec)>(sqrt(M)-1)));
    Ydec(imag(Ydec)<0) = real(Ydec(imag(Ydec)<0));
    Ydec(:) = Ydec - (sqrt(M)/2-0.5)*(1+1j);
    Ydec(pilot_carr) = pilot_data; 
    Ydec(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 0;
    
    % 8) Signal reconstruction reconstrucción del símbolo OFDM en el tiempo
    y_rx(:) = ifft(Ydec);
    % 9) Add cyclic prefix
    y_rx_cp(:) = [y_rx(end-CP+1:end); y_rx];
    % 10) Save reconstructed signal Guardar el símbolo reconstruido
    y_dp(sof(i)+1:sof(i)+CP+T_symb) = y_rx_cp;
    % 11) Obtain Clutter-free observation channel signal
    
    Yobs(:) = Y-Ydec.*Haux;                                         % multiplicacion: recreacion inversa
    y_obs(:) = ifft(Yobs);
    y_obs_cp(:) = [y_obs(end-CP+1:end); y_obs];
    y_oc(sof(i)+1:sof(i)+CP+T_symb) = y_obs_cp;
end

y_oc(isnan(y_oc))=0;        % señal eco
y_dp(isnan(y_dp))=0;        % señal referencia

% y_oc(y_oc==0) = [];
% y_dp(y_dp==0) = [];

figure;
plot(10*log10(abs(y_dp)));
title('Señal referencia');

figure;
plot(10*log10(abs(y_oc)));
title('Señal eco');

%% Añadir ceros en frecuencia

% ventana = rectwin(length(y_oc));
% 
% Y_OC = fftshift(fft(y_oc));
% Y_DP = fftshift(fft(y_dp));
% 
% ventana = [pad ; ventana ; pad];
% 
% Y_OC = [pad ; Y_OC ; pad];
% Y_DP = [pad ; Y_DP ; pad];
% 
% Y_OC = Y_OC .* ventana;
% Y_DP = Y_DP .* ventana;
% 
% y_oc = ifft(ifftshift(Y_OC));
% y_dp = ifft(ifftshift(Y_DP));

%% Interpolacion clasica
% 
% F = 1;
% y_oc_interp = zeros(1,length(y_oc)*F);
% y_dp_interp = zeros(1,length(y_dp)*F);
% 
% y_oc_interp(1:F:end) = y_oc;
% y_dp_interp(1:F:end) = y_dp;
% 
% y_oc = y_oc_interp';
% y_dp = y_dp_interp';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Correlation processing (Batch Algorithm)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% TODO: Graficar posibles valores de L y su efectividad

L = 2^12; % Decimation factor (Batch block processing length)

s_rx = y_oc;
s_tx = y_dp(1:min(length(y_dp),length(s_rx)));
s_tx = y_dp;

%Funcion mia
factor_interp = 10;
[Y2,num_exec] = batch_corr_algorithm(s_rx,s_tx,L,factor_interp);

% Manejo de señales para evitar problemas por longitud

s_rx = [s_rx(:); zeros(L-mod(length(s_rx),L),1)];
s_tx = [s_tx(:); zeros(length(s_rx)-length(s_tx),1)];

% Colocación en grupos de longitud L.

s_rx = reshape(s_rx,L,length(s_rx)/L);
s_tx = reshape(s_tx,L,length(s_tx)/L);

% Ahora las señales son de la misma longitud


% Numero de columnas igual al numero de grupos

cols = size(s_rx,2); 

% Preallocation

Y = zeros(L,cols);

Lfft_aux = 2^(nextpow2(2*L-1));

y_aux = zeros(Lfft_aux,1);
factor_interp = 10;
pad = zeros(((length(s_rx)*(factor_interp -1))/2),1);

window_prev = rectwin(Lfft_aux);
window_prev = [pad ; window_prev ; pad];
    
% Time-domain processing
for i = 1:cols
    y_aux(:) = ifft(fft(s_rx(:,i),Lfft_aux).*conj(fft(s_tx(:,i),Lfft_aux)),Lfft_aux);
    Y(:,i) = y_aux(1:L);
end

% Optimo cuando la longitud es potencia de base 2 (FFT)
K = 2^(nextpow2(cols)); % FFT length
Y = [Y, zeros(L,K-cols)]; % zero-padding

% Funcion de enventanado, reducción de lobulos deterministas
% TODO: Graficar efectos de diferentes funciones de enventanado
% Window function applied to reduce the sidelobe level

% w = chebwin(cols,80);
% w = w(:)/sum(w);
% w = [w; zeros(K-cols,1)];
% 
% Y = Y.*repmat(w.',L,1);
    
% PREGUNTAR FFTSHIFT
Y = fftshift(fft(Y,K,2),2);
    
figure('Name','Alejandro');
a1 = axes();
surf(a1,((0:(2^nextpow2(num_exec)-1))/2^nextpow2(num_exec)-0.5)*(fs/L)*3e8/(2*f),(0:(L-1))*3e8/(2*fs),10*log10((abs(Y2)).^2))
xlim(a1,[-3e8*fs/(4*L*f),3e8*fs/(4*L*f)])
ylim(a1,[0,(L-1)*3e8/(2*fs)])
ylabel(a1,'Bistatic Range (m)','Interpreter','Latex')
xlabel(a1,'Bistatic Velocity (m/s)','Interpreter','Latex')
zlabel(a1,'Power Level (dB)','Interpreter','Latex')
title(a1,['DNR = ',num2str(SNR),' dB, $T_{coh} = ',num2str(num_exec*L/fs,'%4.2f'),'$ sec.'],'Interpreter','Latex')
shading interp
colormap jet

f1 = figure('Name','Carlos');
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


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% SNR
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Windows
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% figure;
% hold on;
% for i=1:5:41
%     plot(10*log10(kaiser(cols,i)),'DisplayName',['kaiser alfa = ' num2str(i)]);
% end
% legend(gca,'show','Location','southeast');
% xlabel('Samples');
% ylabel('Power Level (dB)');
% 
% figure;
% w = hanning(cols);
% plot(linspace(-fs*0.5,fs*0.5,length(w)),(abs(fftshift(fft(w)))));
