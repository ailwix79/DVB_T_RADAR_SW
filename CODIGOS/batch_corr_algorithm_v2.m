
% PASSIVE RADAR TEST

%% SCENARIO GENERATION

SNR = 60;
[s_rx,pilot_carr,pilot_data,M,fs,T_symb,CP,n_symb,f] = scenario_generator(SNR);


%% CHANNEL ESTIMATION AND SIGNAL SEPARATION

s_x = conj([s_rx; zeros(T_symb,1)]).*[zeros(T_symb,1); s_rx];      
s_x = conv(s_x, ones(CP,1)); 

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
    Haux_coax = Haux;
    
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

s_rx = y_oc;
s_tx = y_dp(1:min(length(y_dp),length(s_rx)));
s_tx = y_dp;

% % Algoritmo de batch
L = 2^12;
lfft = 2^(nextpow2(2*L-1));
factor_interp = 4;
    
[Y2,cols] = batch_corr_algorithm(s_rx,s_tx,L,factor_interp);

figure('Name','Alejandro');
a1 = axes();
surf(a1,((0:(2^nextpow2(cols)-1))/2^nextpow2(cols)-0.5)*(fs/L)*3e8/(2*f),(0:(L-1))*3e8/(2*fs),10*log10((abs(Y2)).^2))
xlim(a1,[-3e8*fs/(4*L*f),3e8*fs/(4*L*f)])
ylim(a1,[0,(L-1)*3e8/(2*fs)])
ylabel(a1,'Bistatic Range (m)','Interpreter','Latex')
xlabel(a1,'Bistatic Velocity (m/s)','Interpreter','Latex')
zlabel(a1,'Power Level (dB)','Interpreter','Latex')
title(a1,['SNR = ',num2str(SNR),' dB, $T_{coh} = ',num2str(cols*L/fs,'%4.2f'),'$ sec.'],'Interpreter','Latex')
shading interp
colormap jet

% DVB-T SCENARIO GENERATOR
SNR = 60;
AVG = 4;

[s_rx,~,~,M,fs,Tu,CP,n_symb,f] = scenario_generator_v2(SNR);

% CONTINUAL AND SCATTER PILOT GENERATION
Kmin = 0;
Kmax = 6816;

l = 0:67;
p = 0:floor(Kmax/12);

wk = prbs_dvbt(Kmax+1);

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
                    
cont_pilot_inf = (4/3)*2*(0.5-wk(cont_pilot_carriers+1));

cont_pilot_carriers_aux = cont_pilot_carriers - Kmax/2;
cont_pilot_carriers_aux(cont_pilot_carriers_aux <0) = cont_pilot_carriers_aux(cont_pilot_carriers_aux<0) + Tu;

scatt_pilot_carriers = [12*p(1:end-1);3+12*p(1:end-1);6+12*p(1:end-1);9+12*p(1:end-1)];
scatt_pilot_inf = (4/3)*2*(0.5-wk(scatt_pilot_carriers+1));

scatt_pilot_carriers_aux = scatt_pilot_carriers - Kmax/2;
scatt_pilot_carriers_aux(scatt_pilot_carriers_aux<0) = scatt_pilot_carriers_aux(scatt_pilot_carriers_aux<0) + Tu;

% CHANNEL ESTIMATION
% Perfect synchronization (signal starts at sample zero)
locs = 0:(Tu+CP):length(s_rx);

% Memories initialization
symb_memory = zeros(Tu,AVG);
H_memory = zeros(Tu,AVG);
y_dp = zeros(Tu+CP,n_symb);
y_oc = zeros(Tu+CP,n_symb);   
memory_init_offset = AVG-1;
n = 1;
for i=1:(n_symb+memory_init_offset)
    H = complex(zeros(1,Tu));
    Y = zeros(Tu,1);
    if (i < AVG)
        symb_memory(:,2:end) = symb_memory(:,1:end-1); 
        symb_memory(:,1) = s_rx(locs(i)+CP+1:min(locs(i+1),length(s_rx)));
        X = fft(symb_memory(:,1),Tu);
        scatter_values = zeros(1,size(scatt_pilot_carriers,1));

        for j=1:size(scatt_pilot_carriers,1)
            scatter_values(j) = abs(sum(X(scatt_pilot_carriers_aux(j,:)+1).*scatt_pilot_inf(j,:)'));
        end

        [~,subframe_number] = max(scatter_values);
        Y(cont_pilot_carriers_aux+1) = cont_pilot_inf;
        Y(scatt_pilot_carriers_aux(subframe_number,:)+1) = scatt_pilot_inf(subframe_number,:);

        H(cont_pilot_carriers_aux+1) = X(cont_pilot_carriers_aux+1)./Y(cont_pilot_carriers_aux+1);
        H(scatt_pilot_carriers_aux(subframe_number,:)+1) = X(scatt_pilot_carriers_aux(subframe_number,:)+1)./Y(scatt_pilot_carriers_aux(subframe_number,:)+1);

        H = fftshift(H);
        carriers = find(abs(H)~=0);
        H = interp1(carriers,H(carriers),1:Tu,'linear','extrap');
        H = H(min(carriers):max(carriers));
        H = [zeros(1,round(((Tu-length(H))/2)-1)), H];
        H = [H, zeros(1,Tu-length(H))];
        H(H==0)=1;
        % Channel memory displacement and introduction of new estimation
        H_memory(:,2:end) = H_memory(:,1:end-1);
        H_memory(:,1) = H';
        
    % Once filter is 3/4 full, this "else" conditional will start working
    else
        if i <= n_symb        
            symb_memory(:,2:end) = symb_memory(:,1:end-1); 
            symb_memory(:,1) = s_rx(locs(i)+CP+1:min(locs(i+1),length(s_rx)));
            X = fft(symb_memory(:,1),Tu);
            scatter_values = zeros(1,size(scatt_pilot_carriers,1));

            for j=1:size(scatt_pilot_carriers,1)
                scatter_values(j) = abs(sum(X(scatt_pilot_carriers_aux(j,:)+1).*scatt_pilot_inf(j,:)'));
            end

            [~,subframe_number] = max(scatter_values);
            Y(cont_pilot_carriers_aux+1) = cont_pilot_inf;
            Y(scatt_pilot_carriers_aux(subframe_number,:)+1) = scatt_pilot_inf(subframe_number,:);

            H(cont_pilot_carriers_aux+1) = X(cont_pilot_carriers_aux+1)./Y(cont_pilot_carriers_aux+1);
            H(scatt_pilot_carriers_aux(subframe_number,:)+1) = X(scatt_pilot_carriers_aux(subframe_number,:)+1)./Y(scatt_pilot_carriers_aux(subframe_number,:)+1);

            H = fftshift(H);
            carriers = find(abs(H)~=0);
            H = interp1(carriers,H(carriers),1:Tu,'linear','extrap');
            H = H(min(carriers):max(carriers));
            H = [zeros(1,round(((Tu-length(H))/2)-1)), H];
            H = [H, zeros(1,Tu-length(H))];
            H(H==0)=1;
            % Channel memory displacement and introduction of new estimation
            H_memory(:,2:end) = H_memory(:,1:end-1);
            H_memory(:,1) = H';
            
            Havg = sum(H_memory,2)./AVG;
        else
            % If out ofd symbol, introduce zeros to move the symbols
            % These zeros are not taken into account for channel estimation
            symb_memory(:,2:end) = symb_memory(:,1:end-1);
            symb_memory(:,1) = zeros(Tu,1);
        end
        % FFT of symbol
        X = fft(symb_memory(:,end), Tu);
        % Equalization
        Yeq = X./Havg;
        % QAM decoder
        Ydec = round(Yeq + (sqrt(M)/2-0.5)*(1+1j));
        Ydec(real(Ydec)>(sqrt(M)-1)) = (sqrt(M)-1) + 1j*imag(Ydec(real(Ydec)>(sqrt(M)-1)));
        Ydec(real(Ydec)<0) = 1j*imag(Ydec(real(Ydec)<0));
        Ydec(imag(Ydec)>(sqrt(M)-1)) = 1j*(sqrt(M)-1) + real(Ydec(imag(Ydec)>(sqrt(M)-1)));
        Ydec(imag(Ydec)<0) = real(Ydec(imag(Ydec)<0));
        Ydec(:) = Ydec - (sqrt(M)/2-0.5)*(1+1j);
        Ydec(cont_pilot_carriers+1) = ((4/3)*2*(0.5-cont_pilot_inf)); 
        Ydec(Tu/2-round(Tu*0.05):Tu/2+round(Tu*0.05),:) = 0;
        % Reference signal reconstruction
        y_rx = ifft(Ydec);
        y_rx_cp = [y_rx(end-CP+1:end); y_rx];
        y_dp(:,n) = y_rx_cp;
        % Echo signal reconstruction
        Yobs = X-Ydec.*Havg;
        y_obs = ifft(Yobs);
        y_obs_cp = [y_obs(end-CP+1:end); y_obs];
        y_oc(:,n) = y_obs_cp;
        % Move n index to allocate reconstructed symbols
        n = n+1;
    end
end

y_dp = reshape(y_dp,[],1);
y_oc = reshape(y_oc,[],1);

s_rx = y_oc;
s_tx = y_dp;

% BATCH ALGORITHM
L = 2^12;
lfft = 2^(nextpow2(2*L-1));
factor_interp = 4;
    
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
title(a1,['SNR = ',num2str(SNR),' dB, $T_{coh} = ',num2str(cols*L/fs,'%4.2f'),'$ sec.'],'Interpreter','Latex')
shading interp
colormap jet