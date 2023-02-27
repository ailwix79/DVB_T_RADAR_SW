
% Focal distance
F = 2.5e3; % m

% Transmitter location
r_tx = [-F,0];

% Receiver location
r_rx = [F,0];

% Evaluate isorange curves (ellipses)
x = linspace(-5*F,5*F,1024); y = x;
[X,Y] = meshgrid(x,y);
d_p_tx = sqrt((r_tx(1)-X).^2+(r_tx(2)-Y).^2);
d_p_rx = sqrt((r_rx(1)-X).^2+(r_rx(2)-Y).^2);
R = d_p_tx + d_p_rx;

% Evaluate the gradient of the range (bistatic velocity)
x2 = linspace(-5*F,5*F,16); y2 = x2;
[X2,Y2] = meshgrid(x2,y2);
d_p_tx = sqrt((r_tx(1)-X2).^2+(r_tx(2)-Y2).^2);
d_p_rx = sqrt((r_rx(1)-X2).^2+(r_rx(2)-Y2).^2);
R2 = d_p_tx + d_p_rx;
[DX, DY] = gradient(R2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PLOTS

% Isorange curves
f1 = figure();
a1 = subplot(1,2,1);
hold(a1,'on')
contour(a1,X,Y,R,'showtext','on')
plot(a1,r_tx(1),r_tx(2),'xr','LineWidth',1.2)
text(a1,r_tx(1)+F/5,r_tx(2)+F/5,'Tx')
plot(a1,r_rx(1),r_rx(2),'xr','LineWidth',1.2)
text(a1,r_rx(1)+F/5,r_rx(2)+F/5,'Rx')
xlabel(a1,'x(m)')
ylabel(a1,'y(m)')
grid(a1,'on')
daspect(a1,[1,1,1])
title(a1,'Bistatic Range')
hold(a1,'off')

% Bistatic Range gradient
a2 = subplot(1,2,2);
hold(a2,'on')
contour(a2,X,Y,R,'showtext','on')
plot(a2,r_tx(1),r_tx(2),'xr','LineWidth',1.2)
text(a2,r_tx(1)+F/5,r_tx(2)+F/5,'Tx')
plot(a2,r_rx(1),r_rx(2),'xr','LineWidth',1.2)
text(a2,r_rx(1)+F/5,r_rx(2)+F/5,'Rx')
xlabel(a2,'x(m)')
ylabel(a2,'y(m)')
grid(a2,'on')
daspect(a2,[1,1,1])
quiver(a2,X2,Y2,DX,DY)
title(a2,'Bistatic Range Gradient')
hold(a2,'off')