%% PRACTICA 6 IMPLEMENTACION DE FILTROS EMPLEANDO DFT


%% CARGA DE DATOS
clear all;
load("LPF_5kHz");
[x,fs] = audioread("PDS_P6_3A_LE1_G4.wav");

%% REPRESENTACION DE DATOS

figure;
plot(linspace(-fs*0.5,fs*0.5,length(x)),abs(fftshift(fft(x))));
title('ESPECTRO');
xlabel('Frecuencia (Hz)');

[h,f] = freqz(b,1,50000,fs);
plot(f,20*log(abs(h)));
title("Respuesta en Frecuencia Filtro");
ylabel("Magnitud (dB)");
xlabel("Frecuencia (Hz)");
xlim([0 8000]);


%% DISEÑO DE UN FILTRO FIR
% Tono a eliminar presente en frecuencia 7kHz

y = filter(b,1,x);

figure;
plot(linspace(-fs*0.5,fs*0.5,length(y)),abs(fftshift(fft(y))));
title('ESPECTRO');
xlabel('Frecuencia (Hz)');

sound(y,fs);

%% IMPLEMENTACION USANDO DFT

% muestras de solape. P-1

P = length(b);
L = 500;
N = L-(P-1);

x1 = [x' zeros(1,L-(rem(N,L)))];
n_ejecuciones = length(x)/N;

B = fft(b,L);   % DFT del filtro con longitud L.

for k = 1:n_ejecuciones
    if (k ~= 1)
        % k != 1, se toman las P-1 muestras anteriores hasta el siguiente
        % múltiplo de N.
        m1 = x1((k-1)*N+1-(P-1): k*N);
    else
        % k = 1, primera iteración. Se toma desde 1 a N, añadiendo ceros al
        % principio para hacer el fragmento de longitud L. Esta línea solo
        % se ejecuta en la primera iteración.
        m1 = [zeros(1,P-1), x1((k-1)*N+1 : k*N)];   
    end
    m2 = ifft(fft(m1).*B);          % convolución circular.
    y((k-1)*N+1:k*N) = m2(P:end);   % se descartan las primeras P-1 muestras
end

%% Analisis de resultados

% Apartado B
g = filter(b,1,x);
t = linspace(0,(length(g)-1)/fs,length(g));

figure;
plot(t,y-g);

figure;
hold on;
stem(t,g);
stem(t,y,"*");
xlabel("Tiempo (s)");
title("Señales Y G en el tiempo");
xlim([0.01 0.0101]);
legend('G','Y');

% Apartado C
e = immse(g,y);

% Apartado D
figure;
hold on;
plot(linspace(-fs*0.5,fs*0.5,length(y)),20*log(abs(fftshift(fft(x)))));
plot(linspace(-fs*0.5,fs*0.5,length(y)),20*log(abs(fftshift(fft(g)))));
plot(linspace(-fs*0.5,fs*0.5,length(y)),20*log(abs(fftshift(fft(y)))));
title('ESPECTRO');
legend('X','G','Y');
xlabel('Frecuencia (Hz)');

% Elimina la energia de las componentes de las bandas eliminadas
