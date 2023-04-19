% Channel estimation and signal separation
% - Alejandro Manuel López Gómez
% - Carlos García de la Cueva

clc;
clear all;

% DVB_T SIGNAL AND CARRIERS DECLARATION

% DVB-T signal parameters (8K signal)

snr_db = 60;
[s_rx,~,~,M,fs,T_symb,CP,n_symb,f] = scenario_generator_v2(snr_db);

avg_symb_num = 4;

% Continual pilots (these are fixed by the standart)
pilot_carriers =      [ 0       48      54      87      141     156     192 ...
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
% therefore all possible k values are

% - kmin = 12p          Sub frame 1 mod(l,4) = 0
% - k = 3 + 12p         Sub frame 2 mod(l,4) = 1
% - k = 6 + 12p         Sub frame 3 mod(l,4) = 2
% - kmax = 9 + 12p      Sub frame 4 mod(l,4) = 3

% Hence, p must be a vector that makes sure k values do not leave the range
% [Kmin, Kmax], that is to say that kmin >= Kmin and kmax <= Kmax
% 12p >= Kmin and 9 + 12p <= Kmax
% therefore p_min = 0 and p_max = (Kmax-9)/12

Kmin = 0;
Kmax = 6816;

l = 0:67;
p = 0:floor(Kmax/12);

scatter_carriers = [12*p(1:end-1);3+12*p(1:end-1);6+12*p(1:end-1);9+12*p(1:end-1)];
scatter_carriers_aux = scatter_carriers - Kmax/2;
scatter_carriers_aux(scatter_carriers_aux<0) = scatter_carriers_aux(scatter_carriers_aux<0) + T_symb;

wk = prbs_dvbt(Kmax+1);

% First symbol starts at 8448 samples (perfect time synchronization)
% locs = (T_symb+CP):(T_symb+CP):length(s_rx);

% CSS SYMBOL START ESTIMATOR

avg_value = 4;
s_x = conj([s_rx; zeros(T_symb,1)]).*[zeros(T_symb,1); s_rx];
s_x = conv(s_x, ones(CP,1));
s_x = s_x(1:avg_value*(T_symb+CP));
A = reshape(s_x,T_symb+CP,length(s_x)/(T_symb+CP));
estimate = sum(abs(A).^2,2);                               
[~,peak_index] = max(estimate);

if (peak_index > (T_symb+CP)/2)
    peak_index = peak_index - (T_symb+CP);
end

locs = peak_index:T_symb+CP:length(s_rx);

% CHANNEL ESTIMATION
% 
% H_avg = zeros(T_symb,avg_symb_num);
% symb_mem = zeros(T_symb,avg_symb_num);
% y_dp = zeros(T_symb+CP,length(s_rx)/(T_symb + CP));
% y_oc = zeros(T_symb+CP,length(s_rx)/(T_symb + CP));     
% memory_init_offset = avg_symb_num-1;
% n = 1;
%         
% for i=1:((length(s_rx)/(T_symb + CP))+memory_init_offset)
%     if (i< avg_symb_num)
%         symb = s_rx(locs(i)+CP+1:min(locs(i+1),length(s_rx)));
%         symb_mem(:,2:end) = symb_mem(:,1:end-1);          
%         symb_mem(:,1) = symb;
%  
%         X = fft(symb_mem(:,1), T_symb);
%             
%         subframe_scatter_values = zeros(1,size(scatter_carriers,2));
%         for j=1:size(scatter_carriers,2)
%             subframe_scatter_values(:,j) = (abs(sum(X(scatter_carriers_aux(:,j)+1).*(((4/3)*2*(0.5-wk(scatter_carriers(:,j)+1)))'))))^2;
%         end
%     
%         [~,subframe_number] = max(subframe_scatter_values);
%             
%         carriers = unique(sort([pilot_carriers' ; scatter_carriers(:,subframe_number)]));
%         H = X(carriers+1)./(((4/3)*2*(0.5-wk(carriers+1)))');
%             
%         H_avg(:,2:end) = H_avg(:,1:end-1);
%         H_avg(:,1) = interp1(carriers,H,1:T_symb,'linear','extrap')';
%         
%     else
%         if i <= length(s_rx)/(T_symb + CP)
% %             symb(1:min(T_symb,length(s_rx)-CP)) = s_rx((T_symb+CP)*(i-1)+CP+1:min((CP+T_symb)*(i),length(s_rx)));
%             symb = s_rx(locs(i)+CP+1:min(locs(i+1),length(s_rx)));
%             symb_mem(:,2:end) = symb_mem(:,1:end-1);          
%             symb_mem(:,1) = symb;
%  
%             X = fft(symb_mem(:,1), T_symb);
%             
%             subframe_scatter_values = zeros(1,size(scatter_carriers,2));
%             for j=1:size(scatter_carriers,2)
%                 subframe_scatter_values(:,j) = (abs(sum(X(scatter_carriers_aux(:,j)+1).*(((4/3)*2*(0.5-wk(scatter_carriers(:,j)+1)))'))))^2;
%             end
%     
%             [~,subframe_number] = max(subframe_scatter_values);
%             
%             carriers = unique(sort([pilot_carriers' ; scatter_carriers(:,subframe_number)]));
%             H = X(carriers+1)./(((4/3)*2*(0.5-wk(carriers+1)))');
%             
%             H_avg(:,2:end) = H_avg(:,1:end-1);
%             H_avg(:,1) = interp1(carriers,H,1:T_symb,'linear','extrap')';
%         else
%             symb = zeros(T_symb,1);
%             symb_mem(:,2:end) = symb_mem(:,1:end-1);
%             symb_mem(:,1) = symb;
%         end
%         
%         Haux = sum(H_avg,2)./avg_symb_num;
%         Yeq = (fft(symb_mem(:,end), T_symb))./Haux;
%         
%         % QAM decoder and symbol reconstruction
%         Ydec = round(Yeq + (sqrt(M)/2-0.5)*(1+1j));
%         Ydec(real(Ydec)>(sqrt(M)-1)) = (sqrt(M)-1) + 1j*imag(Ydec(real(Ydec)>(sqrt(M)-1)));
%         Ydec(real(Ydec)<0) = 1j*imag(Ydec(real(Ydec)<0));
%         Ydec(imag(Ydec)>(sqrt(M)-1)) = 1j*(sqrt(M)-1) + real(Ydec(imag(Ydec)>(sqrt(M)-1)));
%         Ydec(imag(Ydec)<0) = real(Ydec(imag(Ydec)<0));
%         Ydec(:) = Ydec - (sqrt(M)/2-0.5)*(1+1j);
%         Ydec(pilot_carriers+1) = ((4/3)*2*(0.5-wk(pilot_carriers+1))); 
%         Ydec(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 0;
%         
%         y_rx = ifft(Ydec);
%         y_rx_cp = [y_rx(end-CP+1:end); y_rx];
%         y_dp(:,n) = y_rx_cp;
% 
%         Yobs = X-Ydec.*Haux;
%         y_obs = ifft(Yobs);
%         y_obs_cp = [y_obs(end-CP+1:end); y_obs];
%         y_oc(:,n) = y_obs_cp;
%         
%         n = n+1;
%     end
% end
% 
% y_dp = reshape(y_dp,[],1);
% y_oc = reshape(y_oc,[],1);

% CHANNEL ESTIMATION V3

H_avg = zeros(T_symb,avg_symb_num);
symb_mem = zeros(T_symb,avg_symb_num);
y_dp = zeros(T_symb+CP,round(length(s_rx)/(T_symb + CP)));
y_oc = zeros(T_symb+CP,round(length(s_rx)/(T_symb + CP)));     
memory_init_offset = avg_symb_num-1;
n = 1;
        
for i=1:((length(s_rx)/(T_symb + CP))+memory_init_offset)
    if (i< avg_symb_num)
        symb = zeros(T_symb,1);
        symb(1:min(T_symb,length(s_rx)-CP)) = s_rx((T_symb+CP)*(i-1)+CP+1:min((CP+T_symb)*(i),length(s_rx)));
        symb_mem(:,2:end) = symb_mem(:,1:end-1);          
        symb_mem(:,1) = symb;
 
        X = fft(symb_mem(:,1), T_symb)';
            
        subframe_scatter_values = zeros(1,size(scatter_carriers,1));
        for j=1:size(scatter_carriers,1)
            subframe_scatter_values(j) = abs(sum(X(scatter_carriers_aux(j,:)+1).*((4/3)*2*(0.5-wk(scatter_carriers(j,:)+1)))))^2;
        end
    
        [~,subframe_number] = max(subframe_scatter_values);
            
        pilot_inf = wk(pilot_carriers+1);
        Y(pilot_carriers+1) = (4/3)*2*(0.5-pilot_inf);

        scatt_inf = wk(scatter_carriers(subframe_number,:)+1);
        Y(scatter_carriers_aux(subframe_number,:)+1) = (4/3)*2*(0.5-scatt_inf);

        H(pilot_carriers+1) = X(pilot_carriers+1)./Y(pilot_carriers+1);
        H(scatter_carriers_aux(subframe_number,:)+1) = X(scatter_carriers_aux(subframe_number,:)+1)./Y(scatter_carriers_aux(subframe_number,:)+1);
        
        carriers = unique(sort([pilot_carriers , scatter_carriers_aux(subframe_number,:)]));
        H_avg(:,2:end) = H_avg(:,1:end-1);
        H_avg(:,1) = interp1(carriers,H(carriers),1:T_symb,'linear','extrap')';
        
    else
        if i <= length(s_rx)/(T_symb + CP)
            symb = s_rx(locs(i)+CP+1:min(locs(i+1),length(s_rx)));
            symb_mem(:,2:end) = symb_mem(:,1:end-1);          
            symb_mem(:,1) = symb;
 
            X = fft(symb_mem(:,1), T_symb);
            
            subframe_scatter_values = zeros(1,size(scatter_carriers,2));
            for j=1:size(scatter_carriers,2)
                subframe_scatter_values(:,j) = (abs(sum(X(scatter_carriers_aux(:,j)+1).*(((4/3)*2*(0.5-wk(scatter_carriers(:,j)+1)))'))))^2;
            end
    
            [~,subframe_number] = max(subframe_scatter_values);
            
            pilot_inf = wk(pilot_carriers+1);
            Y(pilot_carriers+1) = (4/3)*2*(0.5-pilot_inf);

            scatt_inf = wk(scatter_carriers(subframe_number,:)+1);
            Y(scatter_carriers_aux(subframe_number,:)+1) = (4/3)*2*(0.5-scatt_inf);

            H(pilot_carriers+1) = X(pilot_carriers+1)./Y(pilot_carriers+1);
            H(scatter_carriers_aux(subframe_number,:)+1) = X(scatter_carriers_aux(subframe_number,:)+1)./Y(scatter_carriers_aux(subframe_number,:)+1);

            H(H==0) = [];
            carriers = unique(sort([pilot_carriers , scatter_carriers_aux(subframe_number,:)]));
            H_avg(:,2:end) = H_avg(:,1:end-1);
            H_avg(:,1) = interp1(carriers,H,1:T_symb,'linear','extrap')';
        else
            symb = zeros(T_symb,1);
            symb_mem(:,2:end) = symb_mem(:,1:end-1);
            symb_mem(:,1) = symb;
        end
        
        Haux = sum(H_avg,2)./avg_symb_num;
        Yeq = (fft(symb_mem(:,end), T_symb))./Haux;
        
        % QAM decoder and symbol reconstruction
        Ydec = round(Yeq + (sqrt(M)/2-0.5)*(1+1j));
        Ydec(real(Ydec)>(sqrt(M)-1)) = (sqrt(M)-1) + 1j*imag(Ydec(real(Ydec)>(sqrt(M)-1)));
        Ydec(real(Ydec)<0) = 1j*imag(Ydec(real(Ydec)<0));
        Ydec(imag(Ydec)>(sqrt(M)-1)) = 1j*(sqrt(M)-1) + real(Ydec(imag(Ydec)>(sqrt(M)-1)));
        Ydec(imag(Ydec)<0) = real(Ydec(imag(Ydec)<0));
        Ydec(:) = Ydec - (sqrt(M)/2-0.5)*(1+1j);
        Ydec(pilot_carriers+1) = ((4/3)*2*(0.5-wk(pilot_carriers+1))); 
        Ydec(T_symb/2-round(T_symb*0.05):T_symb/2+round(T_symb*0.05),:) = 0;
        
        y_rx = ifft(Ydec);
        y_rx_cp = [y_rx(end-CP+1:end); y_rx];
        y_dp(:,n) = y_rx_cp;

        Yobs = X-Ydec.*Haux;
        y_obs = ifft(Yobs);
        y_obs_cp = [y_obs(end-CP+1:end); y_obs];
        y_oc(:,n) = y_obs_cp;
        
        n = n+1;
    end
end

y_dp = reshape(y_dp,[],1);
y_oc = reshape(y_oc,[],1);