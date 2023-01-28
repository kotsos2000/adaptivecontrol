clear
clc
close all

% initialize matrices and other values needed
A=[0 1 ; -20 -2];
B=[0;2];
Ct=[1 0];
g=1;
F=-1;
c0=0.125;
% meta apo prakseis
T1=1;
T2=-9.5;
T3=9.375;
Ac=[A+B*T3*Ct B*T1 B*T2 ; g*T3*Ct F+(g*T1) g*T2 ; g*Ct 0 F ];
Bc=[B*c0 ; g*c0 ; 0];
Cc=[ Ct 0 0];
tspan=0:0.1:10 ;

% initialize depending on the question

% question=='a' or 'b' or 'c' or 'd' 

question='c' ; 

% isws gia d na min theloyme sto plot to imitono
% an den to theloyme bgazoyme apo to 2o case to d kai to vazoume sto prwto

%  x = [x1 x2 w1 w2 xm1 xm2]
if question=='a'
    x0=[0.1745 0 0 0 0.1745 0 ];
    % we give this dirac as input so that the system is 0 at t=10
    rfun=@(t) dirac(t-10);
end

if question=='b'
    x0=[0.8727 0 0 0 0.8727 0 ];
    % edw parolo poy o pinakas Ac einai eystathis den exoyme 
    % kali simperifora logw tis mi-grammikothtas kai to sistima apeirizei
    % an prosomoiwsoyme to grammiko systima tha doyme oti eimaste kala
    rfun=@(t) dirac(t-10);
end

if question=='c'
    x0=[0 0 0 0 0 0];
    rfun=@(t) 0.0175*sin(0.5*t);
end

if question=='d'
    x0=[0 0 0 0 0 0];
    rfun=@(t) 0.0873*sin(90*t);
end

% ode to simulate the system and plot after
odefun=@(t,x)system(t,x,rfun,Ac,Bc);
[t, x] = ode23s(odefun, tspan, x0);

figure(1)
if (question=='a' || question=='b') 
    plot(t, x(:, 5), 'y', t, x(:, 6), 'k' , t,x(:,1),'b',t,x(:,2),'g');
    legend({'Model Position', 'Model Velocity','System Position', 'System Velocity'});
end

if (question=='c' || question=='d')
    plot(t, x(:, 5), 'y', t, x(:, 6), 'k' , t,x(:,1),'b',t,x(:,2),'g',t,rfun(t),'r');
    legend({'Model Position', 'Model Velocity','System Position', 'System Velocity','Sine'});
end

title('Model+System')

% den eimai katholou sigouros gia ayto
theta=[1.0000 -9.5000 9.3750 0.1250];
u=theta * [x(:,3)  x(:,4)  x(:,1)  rfun(t)].';
figure(2)
plot(t,u,'r')
legend({'U'})

% edw e=0 , den to kanw plot

function DX = system(t, x, rfun,Ac,Bc)
    % xm[xm1 xm2 w1 w2]
    DXm=Ac*[x(5) ; x(6) ; x(3) ; x(4)]+Bc*rfun(t);
    C = [1 0];
    % u in R -> 1x4 * 4x1
    theta=[1.0000 -9.5000 9.3750 0.1250];
    y=C*x(1:2);
    omega=[x(3); x(4); y; rfun(t)];
    u = theta * omega;
    % w1dot = -2w1 + u
    dw1=-2*x(3)+u;
    dw2=-2*x(4)+y;
    dx1=x(2);
    %gia to b an theloyme na doyme tin symperifora toy grammikoy systhmatos 
    %dx2=2*u-2*x(2)-20*x(1);
    dx2=2*u-2*x(2)-20*sin(x(1));
    DX = [dx1; dx2 ;  dw1; dw2 ; DXm(1:2)];
end