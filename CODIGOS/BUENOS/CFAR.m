
% Cell Averaging CFAR
clc;
close all;

Y = ones(100,100);
x = randi(100);
y = randi(100);
Y(x,y) = 10;
Y = Y + 0.1*randn(100);

N = 24;
G = 4;
P_fa = 0.01;
figure;
imagesc(Y);

Y = padarray(Y,[N+G,N+G],mean(Y,"all"),'both');

[R,D] = size(Y);
alpha = N*(P_fa^(-1/N)-1);
cfar = zeros(R,D);
detections = zeros(R,D);

for i = (N+G+1):(R-N-G)
    for j = (N+G+1):(D-N-G)
        train_avg = (sum(sum(Y(i-N-G:i+N+G, j-N-G:j+N+G))) - 2*sum(Y(i,j)) - sum(sum(Y(i-G:i+G, j-G:j+G))))/((2*(N+G))^2-(2*G)^2);
        threshold = alpha * train_avg;
        cfar(i,j) = threshold;
        
        if Y(i,j) > threshold
            detections(i,j) = 1;
        end
    end
end

detections = detections(N+G+1:R-N-G,N+G+1:D-N-G);
cfar = cfar(N+G+1:R-N-G,N+G+1:D-N-G);

figure;
imagesc(detections);
colorbar;

figure;
surf(cfar);


