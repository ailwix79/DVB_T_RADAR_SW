
% Ganancia antena transmisión (dB)
Gt = 10;

% Ganancia antena recepción (dB)
Gr = 10;

 % Frecuencia (Hz)
f = 20e6;

% Longitud de onda (m)
lambda = freq2wavelen(f);

% Potencia transmitida (W)
Pt = 5;

% RCS(m^2)
rcs = 100;

% Muestras rango eje de abscisas (m)
R = 0:100:100000;

Pr = (Pt*(10^(Gt/10))*(10^(Gr/10))*(lambda^2)*rcs)./(((4*pi)^3)*(R.^4));
    
figure;
hold on;
plot(R/1000,10*log10(1000*Pr));
title(['Received Power vs Range.',' Pt = ',num2str(Pt/1e6),' MW ','Gt = ',num2str(Gt), ' dB ','f = ',num2str(f/1e6), ' MHz']);
xlabel('Range (km)');
ylabel('Received Power (dBm)');
grid;




