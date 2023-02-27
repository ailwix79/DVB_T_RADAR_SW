% Demostración efecto de la interpolación en el dominio de la frecuencia

%% PARTIENDO DE UNA SEÑAL ESCALÓN EN EL DOMINIO DE LA FRECUENCIA

N=30;
X = ones(N,1);

L1=4;
L2=8;
L3=16;
L = [L1 ; L2 ; L3];

xL1 = fftshift(ifft([zeros(((length(X)*(L1-1))/2),1) ; ones(N,1) ; zeros(((length(X)*(L1-1))/2),1)]));
xL2 = fftshift(ifft([zeros(((length(X)*(L2-1))/2),1) ; ones(N,1) ; zeros(((length(X)*(L2-1))/2),1)]));
xL3 = fftshift(ifft([zeros(((length(X)*(L3-1))/2),1) ; ones(N,1) ; zeros(((length(X)*(L3-1))/2),1)]));
x = fftshift(ifft(X));

xL1 = xL1./max(xL1);
xL2 = xL2./max(xL2);
xL3 = xL3./max(xL3);
x = x./max(x);

figure;
hold on;
grid on;
xlim([-(N/2)*max(L1,L2) (N/2)*max(L1,L2)-1]);
ylim([-20 0]);
plot(-(N/2)*L1:(N/2)*L1-1,10*log10(abs(xL1)));
plot(-(N/2)*L2:(N/2)*L2-1,10*log10(abs(xL2)));
plot(-(N/2)*L3:(N/2)*L3-1,10*log10(abs(xL3)));
plot(-(N/2):(N/2)-1,10*log10(abs(x)));
title('Rx[n]');
legend([append('L = ',cellstr(int2str(L))) ; 'Señal No Interpolada']);

% La interpolación en frencuencia muestra los lóbulos secundarios, aunque
% diferentes valores de interpolación no modifican o modifican muy poco en nivel o la relación
% entre en lóbulo primario y secundario

%% ENVENTANADO

N=30;
X = ones(N,1);

L1=16;
L2=8;
L3=16;
L = [L1 ; L2 ; L3];

XL1 = [zeros(((length(X)*(L1-1))/2),1) ; ones(N,1) ; zeros(((length(X)*(L1-1))/2),1)];
XL2 = [zeros(((length(X)*(L2-1))/2),1) ; ones(N,1) ; zeros(((length(X)*(L2-1))/2),1)];
XL3 = [zeros(((length(X)*(L3-1))/2),1) ; ones(N,1) ; zeros(((length(X)*(L3-1))/2),1)];

% XL1 = XL1.*hamming(length(XL1));
XL2 = XL2.*hamming(length(XL2));
XL3 = XL3.*hamming(length(XL3));
X = X.*hamming(length(X));

xL1 = fftshift(ifft(XL1));
xL2 = fftshift(ifft(XL2));
xL3 = fftshift(ifft(XL3));
x = fftshift(ifft(X));

xL1 = xL1./max(xL1);
xL2 = xL2./max(xL2);
xL3 = xL3./max(xL3);
x = x./max(x);

figure;
hold on;
grid on;
xlim([-(N/2)*max(L1,L2) (N/2)*max(L1,L2)-1]);
ylim([-20 1]);
plot(-(N/2)*L1:(N/2)*L1-1,10*log10(abs(xL1)));
plot(-(N/2)*L2:(N/2)*L2-1,10*log10(abs(xL2)));
plot(-(N/2)*L3:(N/2)*L3-1,10*log10(abs(xL3)));
plot(-(N/2):(N/2)-1,10*log10(abs(x)));
title('Rx[n]');
legend([append('L = ',cellstr(int2str(L))) ; 'Señal No Interpolada']);

% La interpolación en frencuencia muestra los lóbulos secundarios, aunque
% diferentes valores de interpolación no modifican o modifican muy poco en nivel o la relación
% entre en lóbulo primario y secundario

