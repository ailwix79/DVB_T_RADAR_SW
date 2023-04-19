%% CLEANUP FROM PREVIOUS RUNS

clear all;
close all;
clc;

%% SCENARIO PREPARAION

SNR = 65;
AVG = 8;
L_inf = 6048;
c = 1/sqrt(60);

[s_rx,~,~,original_data,M,fs,Tu,CP,n_symb,f] = scenario_generator_DVBT(SNR);

%% DECLARATION OF FIXED BY STANDART SIGNAL PARAMETERS

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

%% CHANNEL ESTIMATION

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
    
    %Initialization phase, no symbols from this phase are used
    if (i < AVG)
        %Displace memory and add new symbol
        symb_memory(:,2:end) = symb_memory(:,1:end-1); 
        symb_memory(:,1) = s_rx(locs(i)+CP+1:min(locs(i+1),length(s_rx)));
        X = fft(symb_memory(:,1),Tu);
        % OFDM subframe type detection
        scatter_values = zeros(1,size(scatt_pilot_carriers,1));
        for j=1:size(scatt_pilot_carriers,1)
            scatter_values(j) = abs(sum(X(scatt_pilot_carriers_aux(j,:)+1).*scatt_pilot_inf(j,:)'));
        end
        [~,subframe_number] = max(scatter_values);
        % Initial chaneel estimation (incomplete estimation)
        Y(cont_pilot_carriers_aux+1) = cont_pilot_inf;
        Y(scatt_pilot_carriers_aux(subframe_number,:)+1) = scatt_pilot_inf(subframe_number,:);
        H(cont_pilot_carriers_aux+1) = X(cont_pilot_carriers_aux+1)./Y(cont_pilot_carriers_aux+1);
        H(scatt_pilot_carriers_aux(subframe_number,:)+1) = X(scatt_pilot_carriers_aux(subframe_number,:)+1)./Y(scatt_pilot_carriers_aux(subframe_number,:)+1);
        % Complete channel estimation (fill null gaps)
        % Shift to avoid bin allocation problems
        H = fftshift(H);
        % Find reference points for interpolation
        carriers = find(abs(H)~=0);
        % Interpolate the whole useful symbol length
        H = interp1(carriers,H(carriers),1:Tu,'linear','extrap');
        % Keep valid channel estimation
        H = H(min(carriers):max(carriers));
        % Accomodate length so its equal to symbol length
        H = [zeros(1,round(((Tu-length(H))/2)-1)), H];
        H = [H, zeros(1,Tu-length(H))];
        % Zero Force Equalizer
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
        % Channel memory displacement and introduction of new channel estimation
        H_memory(:,2:end) = H_memory(:,1:end-1);
        H_memory(:,1) = H';
        % Calculate averaged response
        Havg = sum(H_memory,2)./AVG;
        % Equalization
        Yeq = X./Havg;
        % QAM decoder
        Ydec = complex(zeros(Tu,1));
        % Prepare signal structure
        data_carriers = zeros(Tu,1);
        data_carriers(data_carriers_indexes+1) = 1;
        data_carriers(cont_pilot_carriers_aux+1) = 0;
        data_carriers(scatt_pilot_carriers_aux(subframe_number,:)+1) = 0;
        data_carriers(tps_carriers+1) = 0;
        % Denormalise QAM values
        data_demmod = round(Yeq(logical(data_carriers))./c);
        % QAM Hard decision decoder
        data_demmod(real(data_demmod)>(sqrt(M))) = (sqrt(M)) + 1j*imag(data_demmod(real(data_demmod)>(sqrt(M))));
        data_demmod(real(data_demmod)<(-1*sqrt(M))) = (-1*sqrt(M)) + 1j*imag(data_demmod(real(data_demmod)<(-1*sqrt(M))));
        data_demmod(imag(data_demmod)>(sqrt(M))) = 1j*(sqrt(M)) + real(data_demmod(imag(data_demmod)>(sqrt(M))));
        data_demmod(imag(data_demmod)<(-1*sqrt(M))) = (-1j*sqrt(M)) + real(data_demmod(imag(data_demmod)<(-1*sqrt(M))));
        % Normalise values
        Ydec(logical(data_carriers)) = data_demmod*c;
        % Add continual and scatter pilots
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

s_rx = y_oc;
s_tx = y_dp;

L = 2^12;
lfft = 2^(nextpow2(2*L-1));
factor_interp = 4;
    
[Y2,cols] = batch_corr_algorithm(s_rx,s_tx,L,factor_interp);
Y2 = Y2-mean(abs(Y2),"all");

figure('Name','Range-Doppler Map');
a1 = axes();
imagesc(a1,((0:(2^nextpow2(cols)-1))/2^nextpow2(cols)-0.5)*(fs/L)*3e8/(2*f),(0:(L-1))*3e8/(2*fs),10*log10((abs(Y2)).^2))
xlim(a1,[-3e8*fs/(4*L*f),3e8*fs/(4*L*f)])
ylim(a1,[0,(L-1)*3e8/(2*fs)])
ylabel(a1,'Bistatic Range (m)','Interpreter','Latex')
xlabel(a1,'Bistatic Velocity (m/s)','Interpreter','Latex')
zlabel(a1,'Power Level (dB)','Interpreter','Latex')
title(a1,['SNR = ',num2str(SNR),' dB, $T_{coh} = ',num2str(cols*L/fs,'%4.2f'),'$ sec.'],'Interpreter','Latex')
shading interp
colormap jet
colorbar

%% CA CFAR ALGORITHM

% Initialization parameters
% CAF values (abs for detection purposes)
Y = abs(Y2);
% Reference cells
N = 20;
% Guard cells
G = 10;
% False alarm probability
P_fa = 0.0001;

% Padding for border detections
Y = padarray(Y,[N+G,N+G],mean(Y,"all"),'both');
[R,D] = size(Y);

% Alpha calculation (for threshold estimation)
alpha = N*(P_fa^(-1/N)-1);
% Final vector initialization
cfar = zeros(R,D);
detections = zeros(R,D);
detected_values = zeros(R,D);

% Initialize loop and limit its span to only valid detections
for i = (N+G+1):(R-N-G)
    if floor(100*i/((R-N-G)-(N+G+1))) > floor(100*(i-1)/((R-N-G)-(N+G+1)))
        clc;
        fprintf('Aplicando CA CFAR ... (%3.0f %%)\n',floor(100*i/((R-N-G)-(N+G+1))));
    end
        
    for j = (N+G+1):(D-N-G)
        % Calculate train average value from reference cells
        train_avg = (sum(sum(Y(i-N-G:i+N+G, j-N-G:j+N+G))) - 2*sum(Y(i,j)) - sum(sum(Y(i-G:i+G, j-G:j+G))))/((2*(N+G))^2-(2*G)^2);
        % Calculate threshold
        threshold = alpha * train_avg;
        % Store threshold value for future plotting
        cfar(i,j) = threshold;
        % Compare detection threshold to CUT
        if Y(i,j) > threshold
            detections(i,j) = 1;
            detected_values(i,j) = Y(i,j);
        end
    end
end

%% BISTATIC PARAMETER ESTIMATION

% Obtain detection groups
[groups,num_groups] = bwlabel(detections);

for i=1:num_groups
    % Isolate detection group
    [r,c] = find(groups==i);
    % Prepare cell indexes
    range = sort(unique(r));
    velocity = sort(unique(c));
    % Extract detection values
    one_det = detected_values(range,velocity);
    % Fit parabola to values
    [pr,Sr] = polyfit(range,sum(one_det,2),2);
    [pv,Sv] = polyfit(velocity,sum(one_det),2);
    % Obtain detection from derivative of polinomial
    iR = roots(polyder(pr));
    iV = roots(polyder(pv));
    
    figure;
    hold on;
    grid on;
    plot(linspace(min(range),max(range)),20*log10(polyval(pr,linspace(min(range),max(range)))));
    plot(range,20*log10(sum(one_det,2)),'-o');
    plot(iR,max(20*log10(polyval(pr,linspace(min(range),max(range))))),'^k','LineWidth',1);
    title("Parabola fitted (Range dimension)");
    xlabel("Range cells");
    ylabel("Amplitude (dB)");
    legend("Fitted parabola","CAF values","Estimated position",'Location','southwest');
    
    figure;
    hold on;
    grid on;
    plot(linspace(min(velocity),max(velocity)),20*log10(polyval(pv,linspace(min(velocity),max(velocity)))));
    plot(velocity,20*log10(sum(one_det)),'-o');
    plot(iV,max(20*log10(polyval(pv,linspace(min(velocity),max(velocity))))),'^k','LineWidth',1);
    title("Parabola fitted (Velocity dimension)");
    xlabel("Velocity cells");
    ylabel("Amplitude (dB)");
    legend("Fitted parabola","CAF values","Estimated position",'Location','southwest');
end

% Eliminate added padding
detections = detections(N+G+1:R-N-G,N+G+1:D-N-G);
cfar = cfar(N+G+1:R-N-G,N+G+1:D-N-G);

figure('Name','Detections');
a1 = axes();
imagesc(a1,((0:(2^nextpow2(cols)-1))/2^nextpow2(cols)-0.5)*(fs/L)*3e8/(2*f),(0:(L-1))*3e8/(2*fs),detections)
xlim(a1,[-3e8*fs/(4*L*f),3e8*fs/(4*L*f)])
ylim(a1,[0,(L-1)*3e8/(2*fs)])
ylabel(a1,'Bistatic Range (m)','Interpreter','Latex')
xlabel(a1,'Bistatic Velocity (m/s)','Interpreter','Latex')
zlabel(a1,'Power Level (dB)','Interpreter','Latex')
title(a1,['SNR = ',num2str(SNR),' dB, $P_{fa} = ',num2str(P_fa),'$'],'Interpreter','Latex')
shading interp
colormap jet
colorbar;

figure('Name','CFAR Thresholds');
a1 = axes();
imagesc(a1,((0:(2^nextpow2(cols)-1))/2^nextpow2(cols)-0.5)*(fs/L)*3e8/(2*f),(0:(L-1))*3e8/(2*fs),cfar)
xlim(a1,[-3e8*fs/(4*L*f),3e8*fs/(4*L*f)])
ylim(a1,[0,(L-1)*3e8/(2*fs)])
ylabel(a1,'Bistatic Range (m)','Interpreter','Latex')
xlabel(a1,'Bistatic Velocity (m/s)','Interpreter','Latex')
zlabel(a1,'Power Level (dB)','Interpreter','Latex')
title(a1,['SNR = ',num2str(SNR),' dB, $P_{fa} = ',num2str(P_fa),'$'],'Interpreter','Latex')
shading interp
colormap jet
colorbar;
