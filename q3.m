% gia s-tropopoiisi dinw tin timi pou thelw sta s1,s2

clear 
clc
close all

% paromoia me prin
A=[0 1 ; -20 -2];
B=[0;2];
Ct=[1 0];
Am=[0 1 ; -0.25 -1];
Bm=[0 ; 0.25];
Cm=[1 0];
tspan=0:0.1:10;

% meta apo prakseis
% K=[-9.875 -0.5] L=[1/8]
% gia na doyme an oi ektimiseis mas einai kales
% se kathe erwtima exw diaforetikes arxikopoihseis twn parametrwn

question='a';

if question=='a'
    x0=[0.1745 0 0.1745 0 0 0 2];
    rfun=@(t) dirac(t-10);
end

if question=='b'
    x0=[0.8727 0 0.8727 0 9 1 10];
    rfun=@(t) dirac(t-10);
end

if question=='c'
    x0=[0 0 0 0 5 2 -2];
    rfun=@(t) 0.0175*sin(0.5*t);
end

if question=='d'
    % edw logw toy grigorou imitonoy arxikopoiw sxetika konta
    % stis alithines alliws to apotelesma tha itan poly kako
    x0=[0 0 0 0 -8 -2 1];
    rfun=@(t) 0.0873*sin(90*t);
end


odefun=@(t,x)system(t,x,rfun,Am,Bm);
[t, x] = ode23s(odefun, tspan, x0);

figure(1)
plot(t, x(:, 3), 'r', t, x(:, 4), 'k' , t,x(:,1),'b',t,x(:,2),'g');
legend({'Model Position', 'Model Velocity','System Position', 'System Velocity'});
title('Model+System')

figure(2)
plot(t,x(:,5),'r',t,x(:,6),'k',t,x(:,7),'g')
legend({'K1','K2','L'})
title('Estimations')

figure(3)
plot(t,x(:,1)-x(:,3),'r',t,x(:,2)-x(:,4),'b')
legend({'Position Error','Velocity Error'})
title('Errors')

function DX = system(t, x, rfun,Am,Bm)
    s1=0;
    s2=0;
    P=1000*[20000000 5000000 ; 5000000 30000000];
    xs=[x(1) ; x(2)];
    xm=[x(3) ; x(4)];
    e=xs-xm;
    DXm=Am*xm + Bm*rfun(t);
    K=[x(5) x(6)];
    L=x(7);
    u=-K*xs+L*rfun(t);
    % - s1*K
    DK=Bm.'*P*e*xs.'- s1*K;
    % - s2*L
    DL=-Bm.'*P*e*rfun(t)- s2*L;
    dx1=x(2);
    dx2=2*u-2*x(2)-20*sin(x(1));
    DX = [dx1; dx2 ;  DXm ; DK(1) ; DK(2) ; DL];
end