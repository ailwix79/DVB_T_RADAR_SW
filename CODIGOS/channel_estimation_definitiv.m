clear all;
clc;

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

% Perfect synchronization (signal starts at sample zero)
locs = 0:(Tu+CP):length(s_rx);

symb_memory = zeros(Tu,AVG);
i = 1;
    H = complex(zeros(Tu,1));
    Y = zeros(Tu,1);
    symb = s_rx(locs(i)+CP+1:min(locs(i+1),length(s_rx)));
    X = fft(symb,Tu);
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
    