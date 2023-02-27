%% Signal Separation Algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Cambiar sincronismo fino, sincronismo grueso, estimación de canal

% Cross-Correlation of the received signal
s_x = conj([s_rx; zeros(T_symb,1)]).*[zeros(T_symb,1); s_rx];      
s_x = conv(s_x, ones(CP,1)); 

W = T_symb+CP;
sof = []; % Comienzo ventana
amp = []; % Final ventana

x_buff = s_x(1:min(T_symb+CP,length(s_x)));         % buffer
[~,r_idx] = max(abs(x_buff));                       % indice del mayor valor del buffer
val = x_buff(r_idx);                                % valor máximo del buffer.
sof(end+1) = r_idx;                                 % inicializar inicio ventana
amp(end+1) = val;                                   % inicializar final ventana


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

y_oc(isnan(y_oc))=0;
y_dp(isnan(y_dp))=0;

int_vec = 1:8:length(y_oc);
y_oc(:) = interp1(y_oc,int_vec,'linear','extrap');
