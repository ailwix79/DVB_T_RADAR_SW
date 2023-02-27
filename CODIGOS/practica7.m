%% PRACTICA 7

%% CARGA DE DATOS
load("PDS_P7_3A_LE1_G4.mat");
[d,fs] = audioread("PDS_P7_3A_LE1_G4_d_n.wav");
[r,~] = audioread("PDS_P7_3A_LE1_G4_r_n.wav");
[x,~] = audioread("PDS_P7_3A_LE1_G4_x_n.wav");

%% OBTENCIION MU MÁXIMA

[r_x,lags] = xcorr(x);

corr_indexes = (0:M)-(0:M)';
corr_indexes = corr_indexes(:);

Rx = zeros(length(corr_indexes),1);

for i = 1:length(corr_indexes)
    Rx(i) = r_x(lags==corr_indexes(i));
end

Rx = reshape(Rx,M+1,M+1);

[U, V] = eig(Rx/Rx(1,1));
mu_max = 1/max(max(V));

%% LMS

N = length(x);

w = zeros(1,M+1);           % coeficientes del filtro
mem = zeros(1,M+1);         % memoria del filtro
e = zeros(1,N);             % señal de error
y = zeros(1,N);             % ruido estimado
coeff = zeros(N,M+1);       % matriz para coeficientes

for i=1:N
    % Desplazar memoria del filtro
    mem(2:end) = mem(1:end-1);
    
    % Actualizar memoria del filtro              
    mem(1) = x(i);
    
    % Aplicacion del filtro
    y(i) = sum(w .* mem);
    
    % Calcular error
    e(i) = d(i) - y(i);
    
    % Almacenar coeficientes usados
    coeff(i,:) = w;
    
    % Estimación nuevos coeficientes
    w = w + 2 * mu * e(i) * mem;
end

%% REPRESENTACION
t = linspace(0,(N-1)/fs,N);

figure;
hold on;
plot(t,d);
plot(t,e);
title("AUDIO ORIGINAL VS FILTRADO");
xlabel("TIEMPO (s)");
legend('d','e');

figure;
hold on;
plot(t,x);
plot(t,y);
title("RUIDO ORIGINAL VS ESTIMADO");
xlabel("TIEMPO (s)");
legend('x','y');
xlim([0.1 0.12]);
err = immse(x,y');

figure;
hold on;
plot(linspace(-fs*0.5,fs*0.5,length(d)),(abs(fftshift(fft(d)))),'b');
plot(linspace(-fs*0.5,fs*0.5,length(e)),(abs(fftshift(fft(e)))),'g');
plot(linspace(-fs*0.5,fs*0.5,length(x)),(abs(fftshift(fft(x)))),'black');
plot(linspace(-fs*0.5,fs*0.5,length(y)),(abs(fftshift(fft(y)))),'r');
title('ESPECTRO');
legend('D','E','X','Y');
xlabel('Frecuencia (Hz)');
xlim([-5000 5000]);

figure;
hold on;
plot(linspace(-fs*0.5,fs*0.5,length(e)),(abs(fftshift(fft(e)))),'g');
plot(linspace(-fs*0.5,fs*0.5,length(y)),(abs(fftshift(fft(y)))),'r');
title('SEÑALES E e Y. COMPARACION PARTE ELIMINADA');
legend('E','Y');
xlabel('Frecuencia (Hz)');
xlim([-5000 5000]);

figure;
hold on;
plot(linspace(-fs*0.5,fs*0.5,length(d)),(abs(fftshift(fft(d)))),'b');
plot(linspace(-fs*0.5,fs*0.5,length(x)),(abs(fftshift(fft(x)))),'black');
title('SEÑALES D y X. SEÑAL ORIGINAL Y RUIDO QUE CONTIENE');
legend('D','X');
xlabel('Frecuencia (Hz)');
xlim([-5000 5000]);

%% COMPROBACION CONVERGENCIA DE FILTROS

[r_dx, lags] = xcorr(d,x);

corr_indexes = (0:M)';

Rdx = zeros(length(corr_indexes),1);
for i = 1:length(corr_indexes)
    Rdx(i) = r_dx(lags==corr_indexes(i));
end

WIENER = inv(Rx)*Rdx;

figure;
hold on;
plot(1:length(x),repmat(WIENER,1,length(x)),'--black');
plot(coeff);
title("CONVERGENCIA FILTRO");
legend('COEFF WIENER');
