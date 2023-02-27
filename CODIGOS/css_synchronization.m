
function [peak_index] = css_synchronization(s_tx_w_n, T_symb, CP, avg_value)
    s_x = conj([s_tx_w_n; zeros(T_symb,1)]).*[zeros(T_symb,1); s_tx_w_n];
    s_x = conv(s_x, ones(CP,1));
    s_x = s_x(1:avg_value*(T_symb+CP));
    A = reshape(s_x,T_symb+CP,length(s_x)/(T_symb+CP));
    estimate = sum(abs(A).^2,2);                               

    % Estimation extraction
    [~,peak_index] = max(estimate);
    % Error calculation
    if (peak_index > (T_symb+CP)/2)
        peak_index = peak_index - (T_symb+CP);
    end



end