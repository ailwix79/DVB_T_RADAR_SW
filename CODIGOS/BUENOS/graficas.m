N = 0:10:1000;
P_fa = 0.0001;
alpha = N.*(P_fa.^(-1./N)-1);

figure;
hold on;
plot(N,alpha);
title(['Alpha VS reference cells P_{fa} = ',num2str(P_fa),'.']);
xlabel("Number of reference cells");
ylabel("Alpha value");