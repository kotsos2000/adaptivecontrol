% se ayto to erwtima efoson exoyme anadrasi eksodoy,den exoyme
% pliroforia gia thn taxitita , opote perimenoyme megali apoklisi
% stin taxythta kai mikroteri stin thesi

%episis sta prwta 2 erwthmata i eksodos mhdenizetai poly grigora
% kai oi ektimhseis mas pane sto 0

clear
clc
close all

% initialize matrices and other values needed
A=[0 1 ; -20 -2];
B=[0;2];
Ct=[1 0];
g=1;
F=-1;
% me to xeri , dn ginotan alliws
Ac=[A(1,1) A(1,2) 0 0 ; A(2,1) A(2,2) 0  0 ; 0 0  F 0 ; 1 0 0 F ];
Bc=[B ; g ; 0];
Cc=[ Ct 0 0];
tspan=0:0.1:10;

% to pragmatiko theta
%theta=[1.0000 -9.5000 9.3750 0.1250]
% gia na gnwrizoyme ek twn ysterwn an kaname kali ektimisi

question='a';

% x= [x1 x2 w1 w2 th1 th2 th3 th4 phi1 phi2 phi3 phi4 xm1 xm2]
if question=='a'
    rfun=@(t) dirac(t-10);
    x0=[0.1745 0 0 0 11 -90 95 2 0 0 0 0 0.1745 0 ];
end

if question=='b'
    rfun=@(t) dirac(t-10);
    x0=[0.8727 0 0 0 -5 20 0 5 0 0 0 0 0.8727 0 ];
end

if question=='c'
    x0=[0 0 0 0 2 -12 13 2 0 0 0 0 0 0 ];
    rfun=@(t) 0.0175*sin(0.5*t);
end

if question=='d'
    x0=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
    rfun=@(t) 0.0873*sin(90*t);
end


odefun=@(t,x)system(t,x,rfun,Ac,Bc,Cc);
[t, x] = ode23s(odefun, tspan, x0);

figure(1)

plot(t, x(:, 1), 'k', t, x(:, 2), 'b',t,x(:,13),'g',t,x(:,14),'r');
legend({'System Position', 'System Velocity','Model Position', 'Model Velocity'});
title('System')



figure(2)
plot(t,x(:,2)-x(:,14),'b',t,x(:,1)-x(:,13),'r')
legend({'Velocity Error','Position Error'})

figure(3)
plot(t, x(:,9), 'k', t, x(:,10), 'r', t, x(:,11),'g', t, x(:,12),'b');
legend({'Theta1', 'Theta2','Theta3','Theta4'});
title('Theta Estimations')

function DX = system(t,x,rfun,Ac,Bc,Cc)
    % declare params
    sigma=0;
    p0=5;
    gamma=diag([ 1000000 2000000 3000000 1000000]);
    C = [1 0];
    % model state equations
    xm=[x(13) x(14) x(3) x(4)];
    DXm=Ac*xm.'+Bc*rfun(t);
    % outputs
    ym=Cc*xm.';
    y=C*x(1:2);
    % error
    epsilon=y-ym;
    % theta vector
    theta=[x(5) x(6) x(7) x(8)];
    % phi vector
    phi=[x(9);x(10);x(11);x(12)];
    % omega vector
    omega=[x(3); x(4); y; rfun(t)];
    % phidot
    Dphi=-p0*phi + omega;
    % thetadot
    Dtheta=-epsilon*gamma*phi -sigma*theta.';
    u = theta * omega+ Dtheta.'*phi;
    dw1=-2*x(3)+u;
    dw2=-2*x(4)+y; 
    dx1=x(2);
    dx2=2*u-2*x(2)-20*x(1);
    % final dot vector
    DX = [dx1; dx2 ;  dw1; dw2 ; Dtheta; Dphi ;DXm(1);DXm(2)];
end
