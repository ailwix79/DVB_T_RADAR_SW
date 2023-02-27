%% Practica 5 PDS
%
%   Alejandro Manuel López Gómez
%   Emilio Kenji Hernández Kuramata
%
%% DATOS COMUNMENTE USADOS. EJECUTAR ESTE BLOQUE SIEMPRE

clear;
data = load("PDS_P5_3A_LE1_G4.mat");

% Filtro original
a = load("LPF_ELLIPTIC").Den;
b = load("LPF_ELLIPTIC").Num;

% Coeficientes Cuantificación
B = 16;
D = 8;

% Puntos freqz
n = 50000;

%% Quantizer
q = quantizer('fixed','round','saturate',[B,D]);

%% CUANTIFICACION DE LOS COEFICIENTES DE UN FILTRO

% Filtro Cuantificado PRIMER FILTRO
a_q = bin2num(q,num2bin(q,a));
b_q = bin2num(q,num2bin(q,b));

%% SECCIONES DE SEGUNDO ORDEN

% Filtro dividido en SOS
[sos,g] = tf2sos(b,a);

% SOS Cuantificado (Directamente)
sos_q = bin2num(q,num2bin(q,sos));
sos_m = [sos_q(1:4) sos_q(5:8) sos_q(9:12) sos_q(13:16) sos_q(17:20) sos_q(21:24)];

%% RAICES EN SECCIONES DE SEGUNDO ORDEN

% Raices (por inspección visual, hay 4 filas y 6 columnas)
b_roots = [roots(sos(1,(1:3))) roots(sos(2,(1:3))) roots(sos(3,(1:3))) roots(sos(4,(1:3)))];
a_roots = [roots(sos(1,(4:6))) roots(sos(2,(4:6))) roots(sos(3,(4:6))) roots(sos(4,(4:6)))];

[b_sos,a_sos] = sos2tf(sos_m,g);

% Raices Cuantificadas (REPARAR GANANCIA + ECM)
a_roots_q = bin2num(q,num2bin(q,a_roots)); % Polos
b_roots_q = bin2num(q,num2bin(q,b_roots)); % Raices

% Raices Cuantificadas filtro original


%% FILTROS A UTILIZAR
% Filtro Original
s1.a = a;
s1.b = b;


% Filtro con coefs. cuantificados
s2.a = a_q;
s2.b = b_q;


% Filtro secciones segundo orden
s3.a = a_sos;
s3.b = b_sos;


% Filtro con raices cuantificadas
s4.a = poly(a_roots_q);
s4.b = poly(b_roots_q);

%% ANÁLISIS GENERAL
%% Apartado A
% a = den. = polos
% b = num. = raices/ceros


[s2.ar,s2.a_ord] = reordenarVector(roots(s1.a),roots(s2.a));
[s2.br,s2.b_ord] = reordenarVector(roots(s1.b),roots(s2.b));

[s3.ar,s3.a_ord] = reordenarVector(roots(s1.a),roots(s3.a));
[s3.br,s3.b_ord] = reordenarVector(roots(s1.b),roots(s3.b));

[s4.ar,s4.a_ord] = reordenarVector(roots(s1.a),roots(s4.a));
[s4.br,s4.b_ord] = reordenarVector(roots(s1.b),roots(s4.b));

% ECM POLOS Y CEROS PRIMER FILTRO (Original)
 ecm_a = ECM([s2.ar; s2.br],[s2.a_ord ; s2.b_ord]);

% ECM POLOS Y CEROS SEGUNDO FILTRO
 ecm_b = ECM([s3.ar; s3.br],[s3.a_ord; s3.b_ord]);
 
% ECM POLOS Y CEROS TERCER FILTRO (EL QUE DEBE TENER MENOS ERROR --> COMPROBAR ASOCIACIONES CORRECTAS ENTRE CADA POLO)
 ecm_c = ECM([s4.ar; s4.br],[s4.a_ord; s4.b_ord]);

% Apartado B

[h1,f1] = freqz(s1.b,s1.a,n,data.Fs);
[h2,f2] = freqz(s2.b,s2.a,n,data.Fs);
[h3,f3] = freqz(s3.b,s3.a,n,data.Fs);
[h4,f4] = freqz(s4.b,s4.a,n,data.Fs);

figure;
hold on;
grid on;
plot(f1,20*log10(abs(h1)));
plot(f2,20*log10(abs(h2)));
plot(f3,20*log10(abs(h3)));
plot(f4,20*log10(abs(h4*g)));
xlabel("FRECUENCIA (Hz)");
ylabel("MAGNITUD (dB)");
title("DIFERENCIAS EN GANANCIA. APARTADO B");
legend("original", "a y b q", "sos q", "raices q");


%Apartado C

figure;
hold on;
grid on;
plot(f1,unwrap(angle(h1)));
plot(f2,unwrap(angle(h2)));
plot(f3,unwrap(angle(h3)));
plot(f4,unwrap(angle(h4)));
xlabel("FRECUENCIA (Hz)");
ylabel("FASE (Grados)");
title("DIFERENCIAS EN FASE. APARTADO C");
legend("original", "a y b q", "sos q", "raices q");
xlim([0 5000]);


%Apartado D

figure;
subplot(2,2,1);
zplane(roots(s1.b),roots(s1.a));
title("original");
subplot(2,2,2);
zplane(roots(s2.b),roots(s2.a));
title("a y b q");
subplot(2,2,3);
zplane(roots(s3.b),roots(s3.a));
title("sos q");
subplot(2,2,4);
zplane(roots(s4.b),roots(s4.a));
title("raices q");