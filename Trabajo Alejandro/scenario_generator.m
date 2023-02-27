
function [s_rx,pilot_carr,pilot_data,M,fs,T_symb,CP,n_symb,f] = scenario_generator(snr_db)
    T_symb  = 8192;                 % Duración de un símbolo, longitud de la IFFT, longitud de un simbolo OFDM sin tener en cuenta CP
    CP      = T_symb/32;            % longitud prefijo cíclico, longitud sym OFDM total T_symb + CP
    n_symb  = 128;                  % numero de símbolos OFDM o portadoras piloto

    M       = 16; % Numero de puntos en la constelación
    symbols = (randi(sqrt(M),[T_symb,n_symb])-(sqrt(M)+1)/2) + 1j*(randi(sqrt(M),[T_symb,n_symb])-(sqrt(M)+1)/2); % Constelación QAM

    pilot_data = randi(2,[256,1])-1.5;                                      % Datos a enviar
    pilot_carr_aux = unique(randi(T_symb,[1,512]),'stable');                % Portadoras (piloto y sub)

    % PREGUNTAR
    pilot_carr = pilot_carr_aux(logical((pilot_carr_aux<(T_symb/2-round(T_symb*0.05)))+(pilot_carr_aux>(T_symb/2+round(T_symb*0.05))))); % pilot carriers

    pilot_carr = pilot_carr(1:256); 	% sub carriers
    symbols(pilot_carr,:) = repmat(pilot_data,1,n_symb); %data modulation M-QAM, each carrier (index) for each data element. The resulting constellation is the transmission

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % High frequency symbols are deleted
    symbols(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 0;     % We are only interested in the lowband signal

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
    snr = 10^(snr_db/10);

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
end