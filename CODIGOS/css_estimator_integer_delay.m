
clear all;
clc;

T_symb  = 8192;                 % Duración de un símbolo, longitud de la IFFT, longitud de un simbolo OFDM sin tener en cuenta CP
CP      = T_symb/32;            % longitud prefijo cíclico, longitud sym OFDM total T_symb + CP
n_symb  = 128;                  % numero de símbolos OFDM o portadoras piloto

min_SNR = -30;
max_SNR = 30;
n_simul = 10;
avg_value = 8;
max_delay = 512;

snr_set = min_SNR:1:max_SNR;
E = zeros(length(snr_set),n_simul);
M = 64;

for i=1:length(snr_set)
    if floor(100*i/length(snr_set)) > floor(100*(i-1)/(length(snr_set)))
        clc;
        fprintf('Realizando simulaciones ... (%3.0f %%)\n',floor(100*i/(length(snr_set))));
    end
    
    symbols = (randi(sqrt(M),[T_symb,n_symb])-(sqrt(M)+1)/2) + 1j*(randi(sqrt(M),[T_symb,n_symb])-(sqrt(M)+1)/2); % Constelación QAM
    pilot_data = randi(2,[256,1])-1.5;                                   
    pilot_carr_aux = unique(randi(T_symb,[1,512]),'stable');              
    pilot_carr = pilot_carr_aux(logical((pilot_carr_aux<(T_symb/2-round(T_symb*0.05)))+(pilot_carr_aux>(T_symb/2+round(T_symb*0.05))))); % pilot carriers
    pilot_carr = pilot_carr(1:256);
    symbols(pilot_carr,:) = repmat(pilot_data,1,n_symb);
    symbols(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 0;

    s_tx = ifft(symbols,T_symb,1);
    s_tx = [s_tx(end-CP+1:end,:); s_tx];
    s_tx = s_tx(:);
    s_tx = s_tx/sqrt(mean(abs(s_tx).^2));
            
        for ii=1:n_simul
            
            % Añadir ruido
%             s_tx_w_n =  s_tx + (1/sqrt(10^(snr_set(i)/10)))*(randn(size(s_tx)) + 1j*randn(size(s_tx)))/sqrt(2);
            s_tx_w_n = awgn(s_tx,snr_set(i),'measured');
            
            % Añadir retardo entero
            n_tram = randi([0,max_delay],1);
            s_tx_w_n = [zeros(n_tram,1) ; s_tx_w_n];
            s_tx_w_n = [s_tx_w_n ; zeros(T_symb+CP-mod(length(s_tx_w_n),T_symb+CP),1)]; 
            
            % Algoritmo
            s_x = conj([s_tx_w_n; zeros(T_symb,1)]).*[zeros(T_symb,1); s_tx_w_n];      
            s_x = conv(s_x, ones(CP,1)); 
            
            s_x = [s_x ; zeros(T_symb+CP-mod(length(s_x),T_symb+CP),1)];
            A = reshape(s_x,T_symb+CP,length(s_x)/(T_symb+CP));
            estimate = sum(abs(A),2);
            [~,peak_index] = max(estimate);
            
            % Error
            
            if (peak_index > (T_symb+CP)/2)
                peak_index = peak_index - (T_symb+CP);
            end
            
            error = n_tram - peak_index;
            
            E(i,ii) = error;
        end
end

RMS_E = rms(E,2);

range = min_SNR:1:max_SNR;

figure;
set(gca,'Yscale','log')
plot(range,RMS_E);
title(['RMS Error (integer delay). Number of simulations per SNR value = ',num2str(n_simul)]);
xlabel('SNR');
