
% Correlation Processing Batch Algorithm.
% Author: Alejandro Manuel López Gómez

function [Y,cols] = batch_corr_algorithm(s_rx,s_tx,L,factor_interp)

    s_rx = [s_rx(:); zeros(L-mod(length(s_rx),L),1)];
    s_tx = [s_tx(:); zeros(length(s_rx)-length(s_tx),1)];
    
    s_rx = s_rx - mean(s_rx);
    s_tx = s_tx - mean(s_tx);
    
    s_rx = reshape(s_rx,L,length(s_rx)/L);
    s_tx = reshape(s_tx,L,length(s_tx)/L);

    cols = size(s_rx,2); 

    Y = zeros(L,cols);

    Lfft_aux = 2^(nextpow2(2*L-1));

    y_aux = zeros(Lfft_aux,1);

    for i = 1:cols
        y_aux = ifft(fft(s_rx(:,i),Lfft_aux).*conj(fft(s_tx(:,i),Lfft_aux)),Lfft_aux);
        y_aux = ifft(ifftshift([zeros(ceil((length(y_aux)*(factor_interp -1))/2),1) ; fftshift(fft(y_aux,length(y_aux))).*factor_interp.*chebwin(length(y_aux)) ; zeros(ceil((length(y_aux)*(factor_interp -1))/2),1)]),length(y_aux));
        Y(:,i) = y_aux(1:L);
    end

    K = 2^(nextpow2(cols));
    Y = [Y, zeros(L,K-cols)];
    w = chebwin(cols);
    w = w(:)/sum(w);
    w = [w; zeros(K-cols,1)];

    Y = Y.*repmat(w.',L,1);
    Y = fftshift(fft(Y,K,2),2);
    % Eliminate continuous component
    Y(1:10,:) = 0;
end