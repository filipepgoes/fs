
function Cruise
clear all
global He Ve betae psie roe zf alphaf m g S c CL0 CLa CLdp CLq...
    CLap CD0 k Cm0 Cma Cmdp Cmq Cma Cmap Cyb Cyda Cydd Clb Clp Clr...
    Clda Cldd Cnb Cnp Cnr Cnda Cndd nv nro Ixx Iyy Izz Ixz...
    roi fmaxi Vi Vs Cnvi mv mr Mi ri CLmax qmax sigma Vs Vq D

    g = 9.80665; m = 120000;

    S = 260; c = 6.61;
    Ixx = 5.55e+6; Iyy = 9.72e+6; Izz = 14.51e+6; Ixz = -3.3e+4;

    CLa =  4.982;  CLdp = .435; CL0 = 0; CLq = -.7; CLap = -.3;
    Cm0 = -.025;   Cma = -1.246; Cmdp = -1.46; Cmap = -5*0; Cmq = -15;
    CD0 = 0.0175;  k = 0.06;
    Cyb = -1.5;    Cyda = .05; Cydd = .3;
    Clb = -1.3;    Clp = -13; Clr = 2.9; Clda = -.33; Cldd = .25;
    Cnb =  1.75;   Cnp = -1.5; Cnr = -7.5; Cnda = -.125; Cndd = -1;

    alphaf = 1*pi/180; zf = 1.5;
    fmaxi = 240000; nv = 0; nro=.75;
    Vi=1; Vs=1; % obs: nv = 0
    Cnvi = 1.8e-5; mv = .62; mr = .122; Mi = .6; ri = .7768; roi = 1.225;

    CLmax = 1; qmax = 5.7e4; 

clc

% ponto inicial
sigma = .8;
Xcruz = sqp([5000 100],@cruzeiro,[],@vinculos);
Hi    = Xcruz(1);
Vi    = Xcruz(2);

% dinamica do cuzeiro
[X0,X] = ode45(@dinam,[0 20000],[Vi 0 Hi 0 m]);

figure(1)
subplot(3,2,1),plot(X(:,4)/60,X(:,1)); xlabel('t (min)'); ylabel('V (m/s)');
subplot(3,2,2),plot(X(:,4)/60,X(:,3)/1000); xlabel('t (min)'); ylabel('H (km)');
subplot(3,2,3),plot(X(:,4)/60,X(:,2)*180/pi); xlabel('t (min)'); ylabel('gama (graus)');
subplot(3,2,4),plot(X(:,4)/60,X(:,5)/1000); xlabel('t (min)'); ylabel('m (ton)');
subplot(3,2,5),plot(X(:,4)/60,X0/1000); xlabel('t (min)'); ylabel('x0 (km)');
subplot(3,2,6),plot(X0/1000,X(:,3)/1000); xlabel('x0 (km)'); ylabel('H (km)');
endfunction

function dfdx = dinam(X0,X)
global He Ve betae psie roe zf alphaf m g S c CL0 CLa CLdp CLq...
    CLap CD0 k Cm0 Cma Cmdp Cmq Cma Cmap Cyb Cyda Cydd Clb Clp Clr...
    Clda Cldd Cnb Cnp Cnr Cnda Cndd nv nro Ixx Iyy Izz Ixz...
    roi fmaxi Vi Vs Cnvi mv mr Mi ri CLmax qmax sigma Vs Vq D...
    Hcruz Vcruz

V = X(1); gama = X(2); H = X(3); t = X(4); massa = X(5);

[ro,drodh,Vsom] = atmosfera(H);
M = V/Vsom;

% condiçoes otimas
aux   = m;
m     = massa;
Xcruz = sqp([H V],@cruzeiro,[],@vinculos);%[5000 100]
m     = aux;
Hcruz = Xcruz(1);
Vcruz = Xcruz(2);

% equilibrio
aux    = m;
m      = massa;
xe0    = [m 0 0]; He = Hcruz; Ve = Vcruz;
xe     = fminsearch(@equil,xe0);
m      = aux;
Fe     = xe(1);
pife   = Fe/Fmax(He,Ve);
alphae = xe(2);
dpe    = xe(3);

CL = CL0 + CLa*alphae + CLdp*dpe;
CD = polar(CL,M);
L  = 1/2 * ro * V^2 * S * CL;
D  = 1/2 * ro * V^2 * S * CD;

% consumo de combustivel
mcp = Fe * Cnvi * (M/Mi)^mv * (ro/ri)^mr;

dVdx     = (Fe*cos(alphae+alphaf) - D - massa*g*sin(gama))/massa     /(V*cos(gama));
dgamadx  = (Fe*sin(alphae+alphaf) + L - massa*g*cos(gama))/(massa*V) /(V*cos(gama));
dHdx     = (V*sin(gama))                                             /(V*cos(gama));
dtdx     = 1                                                         /(V*cos(gama));
dmassadx = -mcp                                                      /(V*cos(gama));

dfdx = [ dVdx ; dgamadx ; dHdx ; dtdx ; dmassadx ];
endfunction

function f = cruzeiro(X)
global m g S CLmax qmax sigma Cnvi mv mr Mi ri Vq Vs D

H = X(1); V = X(2);

[ro,drodh,Vsom] = atmosfera(H);

Vs = sqrt(2*m*g/(ro*S*CLmax));
Vq = sqrt(2*qmax/ro);

M = V/Vsom;

CL = 2*m*g/(ro*V^2*S);
CD = polar(CL,M);
D  = 1/2*ro*V^2*S*CD;

F = D;

mcp = F * Cnvi * (M/Mi)^mv * (ro/ri)^mr;
f = 1e6*(sigma*mcp + 1 - sigma)/V;
endfunction

function c=vinculos(X)
global Vs Vq D
H = X(1); 
V = X(2);

c1 = D - Fmax(H,V); % D < Fmax
c2 = Vs - V;        % Vs < V
c3 = V - Vq;        % V < Vq
c = [-c1 -c2 -c3];
endfunction

% polar de arrasto
function CD=polar(CL,M)
global CD0 k
km=k;
if M<=.6
    CD0m=CD0;
elseif M <=1;
    CD0m=0.25375-0.984375*M+1.3125*M^2-0.546875*M^3;
elseif M>1
    CD0m=.035; km=.12*M-.06;
end
CD=CD0m+km*CL^2;
endfunction

% atmosfera
function [ro,drodh,Vsom]=atmosfera(H)
ro1=(H<=11000).*(1.225*(1-2.2558e-5*H).^4.256060537);
T1=(H<=11000).*(288.15-6.5e-3*H);
Vsom1=(H<=11000).*(20.0465*sqrt(T1));
drodh1=(H<=11000).*(-1.176078e-4*(1-2.25577e-5*H)^3.256);
ro2=(H>11000).*(0.3638879*exp(-1.576939464e-4*(H-11000)));
Vsom2=(H>11000).*(295.065141);
drodh2=(H>11000).*(-6.39228e-5*exp(-1.7567e-4*(H-11000)));
ro=ro1+ro2;
Vsom=Vsom1+Vsom2;
drodh=drodh1+drodh2;
endfunction

% traçao maxima do motor
function y=Fmax(H,V)  
global roi fmaxi Vi Vs nv nro
ro=atmosfera(H);
ros=atmosfera(11000);
fmaxs=fmaxi*(ros/roi).^nro.*(Vs/Vi).^nv;
y1=(H<=11000).*(fmaxi*(ro/roi).^nro.*(V/Vi).^nv);
y2=(H>11000).*(fmaxs*(ro/ros).*(V/Vs).^nv);
y=y1+y2;
endfunction


%..........................................................................

function f = equil(xe)
global He Ve zf alphaf m g S c CL0 CLa CLdp CD0 k Cm0 Cma Cmdp

Fe = xe(1);
alphae = xe(2);
dpe = xe(3);

ro = atmosfera(He);

CLe = CL0 + CLa * alphae + CLdp * dpe;
CDe = CD0 + k * CLe^2;
Cme = Cm0 + Cma*alphae + Cmdp*dpe;

De = 1/2 * ro * Ve^2 * S * CDe;
Le = 1/2 * ro * Ve^2 * S * CLe;
Ma = 1/2 * ro * Ve^2 * S * c * Cme;
Mf = Fe * zf;

f1 = Fe*cos(alphae + alphaf) - De;
f2 = Fe*sin(alphae + alphaf) + Le - m*g;
f3 = Ma + Mf;

f = f1^2 + f2^2 + f3^2;
endfunction

Cruise()