% This function provides a pseudorandom binary sequence used to modulate
% continual and scatter pilots in DVB-T

function [wk] = prbs_dvbt(len_seq)

wk = zeros(1,len_seq);

mem = (ones(1,11));

for i = 1:len_seq
    wk(i) = mem(1);
    mem(:) = [mem(2:end) bitxor(mem(1),mem(3))];
end