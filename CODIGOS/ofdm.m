
% GENERACIÓN DE UNA SEÑAL OFDM

M = 16;           % Orden de la modulación QAM
nfft = 128;       % Numero de símbolos OFDM o portadoras
Tsym = 256;       % Subportadoras por simbolo OFDM
cp = 128;      % Longitud prefijo circular

dataIn = randi([0 M-1],nfft,Tsym);
qamSig = qammod(dataIn,M,'UnitAveragePower',true);
y = ofdmmod(qamSig,nfft,cp);

figure;
plot(10*log10(abs(y)));
title('Generated OFDM spectrum');

y_corr = conj([y; zeros(Tsym,1)]).*[zeros(Tsym,1); y];   
y_corr = conv(y_corr, ones(cp,1)); 

figure;
plot(10*log10(abs(ifft(y))));
title('Generated OFDM time domain signal');

figure;
plot(abs(y_corr));
title('Cross correlation');

