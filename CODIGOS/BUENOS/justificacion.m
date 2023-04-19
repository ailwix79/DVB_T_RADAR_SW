% Caracterización algoritmo de Batch
% Demostración de la efectividad de la necesidad de la interpolación y el
% enventanado previos al algoritmo de batch



% EFECTO DE LA INTERPOLACIÓN

SNR = 60;
n = 200;
x = awgn(ones(n,1),SNR);

% Interpolación en frecuencia de la señal
factor_interp = [1,8,16];
X1 = [zeros(ceil((length(x)*(factor_interp(1) -1))/2),1) ; x ; zeros(ceil((length(x)*(factor_interp(1) -1))/2),1)];
X2 = [zeros(ceil((length(x)*(factor_interp(2) -1))/2),1) ; x ; zeros(ceil((length(x)*(factor_interp(2) -1))/2),1)];
X3 = [zeros(ceil((length(x)*(factor_interp(3) -1))/2),1) ; x ; zeros(ceil((length(x)*(factor_interp(3) -1))/2),1)];

x1 = fftshift(ifft(X1));
x2 = fftshift(ifft(X2));
x3 = fftshift(ifft(X3));

x1 = x1./max(x1);
x2 = x2./max(x2);
x3 = x3./max(x3);

figure;
hold on;
grid on;
ylim([-30 0]);
xlim([-100 100]);
plot(-(n/2)*factor_interp(1):(n/2)*factor_interp(1)-1,20*log10(abs(x1)));
plot(-(n/2)*factor_interp(2):(n/2)*factor_interp(2)-1,20*log10(abs(x2)));
plot(-(n/2)*factor_interp(3):(n/2)*factor_interp(3)-1,20*log10(abs(x3)));
legend(append('L = ',cellstr(int2str(factor_interp'))));
title('Rx[n] Original');
xlabel('Muestras');
ylabel('Amplitud (dB)');



% APLICACIÓN VENTANA EN FRENCUENCIA
% Partiendo de una señal AWGN
SNR = 60;
n = 200;
x = awgn(ones(n,1),SNR);

% Interpolación en frecuencia de la señal
factor_interp = [1,8,16];
X1 = [zeros(ceil((length(x)*(factor_interp(1) -1))/2),1) ; x.*chebwin(length(x)) ; zeros(ceil((length(x)*(factor_interp(1) -1))/2),1)];
X2 = [zeros(ceil((length(x)*(factor_interp(2) -1))/2),1) ; x.*chebwin(length(x)) ; zeros(ceil((length(x)*(factor_interp(2) -1))/2),1)];
X3 = [zeros(ceil((length(x)*(factor_interp(3) -1))/2),1) ; x.*chebwin(length(x)) ; zeros(ceil((length(x)*(factor_interp(3) -1))/2),1)];

x1 = fftshift(ifft(X1));
x2 = fftshift(ifft(X2));
x3 = fftshift(ifft(X3));

x1 = x1./max(x1);
x2 = x2./max(x2);
x3 = x3./max(x3);

figure;
hold on;
grid on;
ylim([-90 0])
xlim([-400 400]);
plot(-(n/2)*factor_interp(1):(n/2)*factor_interp(1)-1,20*log10(abs(x1)));
plot(-(n/2)*factor_interp(2):(n/2)*factor_interp(2)-1,20*log10(abs(x2)));
plot(-(n/2)*factor_interp(3):(n/2)*factor_interp(3)-1,20*log10(abs(x3)));
legend(append('L = ',cellstr(int2str(factor_interp'))));
title('Rx[n] Con Enventanado en Frecuencia');
xlabel('Muestras');
ylabel('Amplitud (dB)');

% demostración modulación
t = -5:0.01:5;
n = -5:5;
x = sinc(t);
y = sinc(n);
x = x./max(x);
y = y./max(y);

figure;
hold on;
grid on;
xlim([-6 6]);
plot(t,x);
stem(n,y,'oR');
title('Rx[n] AWGN');
xlabel('Muestras');
ylabel('Amplitud Normalizada');
legend('Rx[n] Verdadera','Rx[n] Muestreada');

z = sinc(t/2);
w = sinc(n/2);

z = z./max(z);
w = w./max(w);

figure;
hold on;
grid on;
xlim([-6 6]);
plot(t,z);
stem(n,w,'oR');
title('Rx[n] AWGN Interpolada (L = 2)');
xlabel('Muestras');
ylabel('Amplitud Normalizada');
legend('Rx[n] Verdadera','Rx[n] Muestreada');

% demostración ventana

n = -30:30;
x = [zeros(20,1) ; ones(21,1) ; zeros(20,1)];
y = fftshift(ifft(x));
y = y./max(y);
x = x./max(x);

h = hamming(61);
z = fftshift(ifft(z,61));

h = h./max(h);
z = z./max(z);
figure;
hold on;
%plot(n,abs(x));
plot(n,abs(y));
% plot(n,abs(h));
plot(n,abs(z));