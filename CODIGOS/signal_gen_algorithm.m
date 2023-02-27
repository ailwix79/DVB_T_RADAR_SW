% OFDM generation algorithm

% Signal Generation

longitud_simbolo = 512;                                                % Longitud simbolo OFDM
datos = randi(2,[256,1]) - 1;                                          % Generacion de datos binarios
banco_portadoras = unique(randi(longitud_simbolo,[1,1024]),'stable');  % Generación de sub carriers y pilotos

% Hay que asignar a cada bit una sub-carrier, las subcarriers sobrantes
% seran portadoras piloto

aux = banco_portadoras(logical((banco_portadoras<(longitud_simbolo/2-round(longitud_simbolo*0.05)))+(banco_portadoras>(longitud_simbolo/2+round(longitud_simbolo*0.05)))));

portadoras = banco_portadoras(1:256);

