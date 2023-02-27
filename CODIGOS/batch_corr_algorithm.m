
% Correlation Processing Batch Algorithm.
% Author: Alejandro Manuel López Gómez
%
%       INPUT:
%       s_rx(n): echo signal (target signal)
%       s_tx(n): reference signal
%       L: Batch block length.
%       factor_interp: Interpolation factor
%
%       OUTPUT:
%       Y: CAF matrix (Cross Ambiguity Function matrix)
%       
%

% function [Y,cols] = batch_corr_algorithm(s_rx,s_tx,L,factor_interp,window_effect)
% 
%     switch window_effect
%         case 'DISTANCE'
%             s_rx = ifft(ifftshift([zeros(ceil((length(s_rx)*(factor_interp -1))/2),1) ; fftshift(fft(s_rx,length(s_rx))).*factor_interp.*chebwin(length(s_rx)) ; zeros(ceil((length(s_rx)*(factor_interp -1))/2),1)]),length(s_rx));
%             s_tx = ifft(ifftshift([zeros(ceil((length(s_tx)*(factor_interp -1))/2),1) ; fftshift(fft(s_tx,length(s_tx))).*factor_interp.*chebwin(length(s_tx)) ; zeros(ceil((length(s_tx)*(factor_interp -1))/2),1)]),length(s_tx));
%         case 'VELOCITY'
%             s_rx = ifft(ifftshift([zeros(ceil((length(s_rx)*(factor_interp -1))/2),1) ; fftshift(fft(s_rx,length(s_rx))) ; zeros(ceil((length(s_rx)*(factor_interp -1))/2),1)]),length(s_rx));
%             s_tx = ifft(ifftshift([zeros(ceil((length(s_tx)*(factor_interp -1))/2),1) ; fftshift(fft(s_tx,length(s_tx))) ; zeros(ceil((length(s_tx)*(factor_interp -1))/2),1)]),length(s_tx));
%             
%             s_rx = s_rx .* factor_interp .* chebwin(length(s_rx));
%             s_tx = s_tx .* factor_interp .* chebwin(length(s_tx));
%         case 'NO'
%     end
% 
%     s_rx = [s_rx(:); zeros(L-mod(length(s_rx),L),1)];
%     s_tx = [s_tx(:); zeros(length(s_rx)-length(s_tx),1)];
%     
%     s_rx = reshape(s_rx,L,length(s_rx)/L);
%     s_tx = reshape(s_tx,L,length(s_tx)/L);
% 
%     cols = size(s_rx,2); 
% 
%     Y = zeros(L,cols);
% 
%     Lfft_aux = 2^(nextpow2(2*L-1));
% 
%     y_aux = zeros(Lfft_aux,1);
% 
%     for i = 1:cols
%         y_aux(:) = ifft(fft(s_rx(:,i),Lfft_aux).*conj(fft(s_tx(:,i),Lfft_aux)),Lfft_aux);
%         Y(:,i) = y_aux(1:L);
%     end
% 
%     K = 2^(nextpow2(cols));
%     Y = [Y, zeros(L,K-cols)];
% 
%     Y = fftshift(fft(Y,K,2),2);
% end


% Correlation Processing Batch Algorithm.
% Author: Alejandro Manuel López Gómez
% 
%       INPUT:
%       s_rx(n): echo signal (target signal)
%       s_tx(n): reference signal
%       L: Batch block length.
%       factor_interp: Interpolation factor
% 
%       OUTPUT:
%       Y: CAF matrix (Cross Ambiguity Function matrix)
      


function [Y,cols] = batch_corr_algorithm(s_rx,s_tx,L,factor_interp,window_effect)

    switch window_effect
        case 'DISTANCE'
            s_rx = ifft(ifftshift(fftshift(fft(s_rx,length(s_rx))).*chebwin(length(s_rx))),length(s_rx));
            s_tx = ifft(ifftshift(fftshift(fft(s_tx,length(s_tx))).*chebwin(length(s_tx))),length(s_rx));
        case 'VELOCITY'
            s_rx = s_rx .* chebwin(length(s_rx));
            s_tx = s_tx .* chebwin(length(s_tx));
        case 'NO'
    end

    s_rx = [s_rx(:); zeros(L-mod(length(s_rx),L),1)];
    s_tx = [s_tx(:); zeros(length(s_rx)-length(s_tx),1)];
    
    s_rx = reshape(s_rx,L,length(s_rx)/L);
    s_tx = reshape(s_tx,L,length(s_tx)/L);

    cols = size(s_rx,2); 

    Y = zeros(L,cols);

    Lfft_aux = 2^(nextpow2(2*L-1));

    y_aux = zeros(Lfft_aux,1);

    for i = 1:cols
        y_aux(:) = ifft(fft(s_rx(:,i),Lfft_aux).*conj(fft(s_tx(:,i),Lfft_aux)),Lfft_aux);
        Y(:,i) = y_aux(1:L);
    end

    K = 2^(nextpow2(cols));
    Y = [Y, zeros(L,K-cols)];

    Y = fftshift(fft(Y,K,2),2);
end