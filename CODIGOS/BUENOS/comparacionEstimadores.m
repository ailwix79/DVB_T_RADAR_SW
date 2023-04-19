% grafica BER OFDM

AVG = 16;
SNR = -20:5:60;
error_OFDM = zeros(length(SNR),1);
error_DVBT = zeros(length(SNR),1);

for i=1:length(SNR)
    [s_tx_original_OFDM, s_tx_reconstruida_OFDM] = graficaBER_OFDM(SNR(i));
    s_tx_original_OFDM(end+1) = 0;
    error_OFDM(i,:) = rms(s_tx_original_OFDM-s_tx_reconstruida_OFDM);
    
    [s_tx_original_DVBT, s_tx_reconstruida_DVBT] = graficaBERDVBT(SNR(i),AVG);
    s_tx_reconstruida_DVBT(end+1) = 0;
    error_DVBT(i,:) = rms(s_tx_original_DVBT-s_tx_reconstruida_DVBT);
end

figure;
subplot(2,1,1);
plot(SNR,error_OFDM);
title("RMSE Previous estimator")
xlabel("SNR (dB)");
ylabel("RMSE");

subplot(2,1,2);
plot(SNR,error_DVBT);
title("RMSE New estimator")
xlabel("SNR (dB)");
ylabel("RMSE");

sgtitle("RMSE reference signal (reconstructed vs original)");

function [s_tx_original_OFDM,s_tx_reconstruida_OFDM] = graficaBER_OFDM(SNR)


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

    s_tx_original_OFDM = s_tx;


    % Direct Path to noise power ratio
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
    
    s_tx_reconstruida_OFDM = y_dp(1:min(length(y_dp),length(s_rx)));

end


function [s_tx_original_DVBT, s_tx_reconstruida_DVBT] = graficaBERDVBT(SNR,AVG)
    % DVB-T SCENARIO GENERATOR
    L_inf = 6048;
    c = 1/sqrt(60);

    [s_rx,s_tx,~,~,M,~,Tu,CP,n_symb,~] = scenario_generator_DVBT(SNR);
    s_tx_original_DVBT = s_tx;
    
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

    tps_carriers = tps_carriers - Kmax/2;
    tps_carriers(tps_carriers<0) = tps_carriers(tps_carriers<0) + Tu;

    data_carriers_indexes = (0:Kmax)-Kmax/2;
    data_carriers_indexes(data_carriers_indexes<0) = data_carriers_indexes(data_carriers_indexes<0) + Tu;

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
    n = 1;

    for i=1:n_symb
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

            % Equalization
            Yeq = X./Havg;

            % QAM decoder
            Ydec = complex(zeros(Tu,1));

            data_carriers = zeros(Tu,1);
            data_carriers(data_carriers_indexes+1) = 1;
            data_carriers(cont_pilot_carriers_aux+1) = 0;
            data_carriers(scatt_pilot_carriers_aux(subframe_number,:)+1) = 0;
            data_carriers(tps_carriers+1) = 0;

            data_demmod = round(Yeq(logical(data_carriers))./c);

            data_demmod(real(data_demmod)>(sqrt(M))) = (sqrt(M)) + 1j*imag(data_demmod(real(data_demmod)>(sqrt(M))));
            data_demmod(real(data_demmod)<(-1*sqrt(M))) = (-1*sqrt(M)) + 1j*imag(data_demmod(real(data_demmod)<(-1*sqrt(M))));
            data_demmod(imag(data_demmod)>(sqrt(M))) = 1j*(sqrt(M)) + real(data_demmod(imag(data_demmod)>(sqrt(M))));
            data_demmod(imag(data_demmod)<(-1*sqrt(M))) = (-1j*sqrt(M)) + real(data_demmod(imag(data_demmod)<(-1*sqrt(M))));

            Ydec(logical(data_carriers)) = data_demmod*c;

            Ydec(cont_pilot_carriers_aux+1) = cont_pilot_inf;
            Ydec(scatt_pilot_carriers_aux(subframe_number,:)+1) = scatt_pilot_inf(subframe_number,:); 

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

    s_rx_reconstructed = y_oc;
    s_tx_reconstruida_DVBT = y_dp;

end