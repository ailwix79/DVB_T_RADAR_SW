
function [eco,ref] = syn_and_recons_algorithm(s_rx,T_symb)
    s_x = conj([s_rx; zeros(T_symb,1)]).*[zeros(T_symb,1); s_rx];      
    s_x = conv(s_x, ones(CP,1)); 

    figure;
    plot(10*log10(abs(s_x)));
    title('Cross correlation');

    % Detect cross correlation peaks (pilot carrier position)
    [~,ind] = max(1:(T_symb + CP));                 % find first pilot carrier index
    I = ind:(T_symb + CP):ind*n_symb; I = I';       % all pilots are separated T_symb + CP
    V = s_rx(I);                                    % get pilot amplitude

end
