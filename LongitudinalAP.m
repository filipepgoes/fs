function LongitudinalAP
clear all;
global g S m c Iyy Cla Cld Cl0 Clq Clap Cm0 Cma Cmd Cmap Cmq Cd0 k...
    nv nro alfaf zf He Ve A B K ictr Xe Ue C ax ah nx nh x0i x0f ivento...
    Kr KI...
    U % In Octave
 
clc;
 
disp('        MOVIMENTO LONGITUDINAL');
disp ' '
 
g = 9.80665; S = 260; m = 120000; c = 6.61; Iyy = 9.72e+6;
 
Cla = 4.982; Cld = .435; Cl0 = 0; Clq = -.7; Clap = -.3; Cm0 = -.025;
Cma = 0; Cmd = -1.46; Cmap = -5; Cmq = -15; Cd0 = 0.0175; k = 0.06;
Cma = -1.246;

nv = 0; nro=.75; alfaf = 1*pi/180; zf = 1.5;
 
% vento
ax = 0; ah = 0;
nx = 1; nh = 1;
x0i = 2000; x0f = 4000;
ivento = 1;
 
% entrada
        He = 7000;
        Ve = 220;
        T = 2*60;
 
disp('Condiçoes de equilibrio:');
He = input (' Altitude de voo (m): ');
Ve = input (' Velocidade de voo (m): ');
disp ' '
        
        dV = 0;
        dgama = 0*pi/180;
        dalfa = 0*pi/180;
        dq = 0;
        dH = 0;
         
        ddp = 0*pi/180;
        dpif = .0;
         
% calculando o equilibrio
xe0 = [m 0 0];
xe = fminsearch(@equil,xe0);
disp('Tracao de equilibrio (N) e posicao da manete');
Fe = xe(1)
pife = Fe/Fmax(He,Ve)
disp('Angulo de ataque e deflexao do profundor de equilibrio (graus)');
xe(2:3)*180/pi
alfae = xe(2);
dpe = xe(3);
 
ictr = 0;
% linearizacao numerica
Xe = [Ve 0 alfae 0 He 0];
Ue = [dpe;pife];
eps1 = 1e-6;
for i=1:5
    X = Xe; X(i)=Xe(i)+eps1; 
    % In Matlab:
    % f1 = dinam(0,X,Ue); 
    % In Octave:
    U = Ue;
    f1 = dinam(0,X); 
    f1(6)=[];
    X = Xe; X(i)=Xe(i)-eps1; 
    % In Matlab:
    % f2 = dinam(0,X,Ue); 
    % In Octave:
    U = Ue;
    f2 = dinam(0,X); 
    f2(6)=[];
    A(:,i) = (f1 - f2) / (2*eps1);
end
 
ilin = input('Deseja ver informaçoes sobre a dinamica linearizada (sim=1)? ')
if ilin == 1
    disp('Linearizaçao');
    disp('Matriz de Estados');
    disp(A);
    damp(A);
    disp('Periodo curto');
    Ac = A(3:4,3:4);
    damp(Ac);
    disp('Raizes');
    raiz = eig(Ac); disp(raiz);
    disp('wn');
    wnc = sqrt(raiz(1)*raiz(2)); disp(wnc);
    disp('Amortecimento');
    ksic=-(raiz(1)+raiz(2))/(2*wnc); disp(ksic);
    roe = atmosfera(He);
    disp('nza');
    nza=1/2*roe*Ve^2*S*Cla/(m*g); disp(nza);
    disp('Periodo longo');
    Al = [A(1:2,1:2) A(1:2,5) ; A(5,1:2) A(5,5)]
    damp(Al);
end
 
for i = 1:2
    U = Ue; U(i) = U(i)+eps1; 
    % In Matlab:
    % f1 = dinam(0,Xe,U); 
    % In Octave:
    f1 = dinam(0,Xe); 
    f1(6)=[];
    U = Ue; U(i) = U(i)-eps1; 
    % In Matlab:
    % f2 = dinam(0,Xe,U); 
    % In Octave:
    f2 = dinam(0,Xe); 
    f2(6)=[];
    B(:,i) = (f1 - f2) / (2*eps1);
end
 
if ilin == 1
    disp('Matriz de Controle');
    B
    Bc = B(3:4,1);
    Bl = [B(1:2,:) ; B(5,:)];
end
 
% Piloto automatico
H = [1 0 0 0 0; 0 0 0 0 1]; F = 0; G = 1;
Achap = [A zeros(5,2); -G*H zeros(2,2)];
Bchap = [B; zeros(2,2)];
 
isim = 1;
while isim ==1
     
    ictr = input('Fazer projeto de um Piloto Automatico por LQR? (1=sim) ');
    K = [];
    if ictr == 1
        % LQR
        Q(6,6)=1; Q(7,7)=1;
        R(1,1)=10; R(2,2)=1e4;
        [K,s,p] = lqr(Achap,Bchap,Q,R);
        disp('Ganhos encontrados:')
        Kr = K(:,1:5)
        KI = -K(:,6:7)
        disp('Polos de malha fechada:')
        disp(eig(Achap-Bchap*K))
    end
     
    disp('Perturbacao na condicao de equilibrio (vento):');
    x0i = input(' Posiçao inicial da perturbaçao (m): ');
    x0f = input(' Posiçao final da perturbaçao (m): ');
    ax  = input(' Velocidade do vento na direçao de voo (m/s): ');
    ah  = input(' Velocidade do vento na direçao ascendete (m/s): ');
    nx  = input(' Numeros de picos na direçao de voo: ');
    nh  = input(' Numeros de picos na direçao ascendente: ');
    ivento = 1;
    disp ' '
%     disp('Perturbacao na condicao de equilibrio em t=0:')
%     dV = input(' dV = ');
%     dgama = input(' dgama (grau)= '); dgama = dgama*pi/180;
%     dalfa = input(' dalfa (grau)= '); dalfa = dalfa*pi/180;
%     dq = input(' dq = ');
%     dH = input(' dH = ');
%     ddp = input(' ddp (grau)= '); ddp = ddp*pi/180;
%     dpif = input(' dpif (%)= '); dpif = dpif/100;
    T = input ('Tempo de simulaçao (mim): '); T = T*60;
    disp ' '
 
    % simulaçao da dinamica longitudinal sem PA
    aux = ictr; ictr = 0;
    X0 = [Ve+dV 0+dgama alfae+dalfa 0+dq He+dH 0];
    U = [dpe+ddp;pife+dpif];
    [t,X] = ode45(@dinam,[0 T],X0,[],U);
    ictr = aux; aux = [];
 
    % simulaçao da dinamica longitudinal com PA ou linear
    if ictr == 1
        X0 = [Ve+dV 0+dgama alfae+dalfa 0+dq He+dH 0 0 0 0 0];
        % In Matlab:
        % [t,X] = ode45(@dinam,[0 T],X0,[],U0);
        % In Octave:
        [t2,X2] = ode45(@dinam,[0 T],X0);
    else
        x0 = [dV dgama dalfa dq dH];
        u = [ddp;dpif];
        % In Matlab:
        % [t2,X2] = ode45(@dinam,[0 T],x0,[],u0);
        % In Octave:
        U = u;
        [t2,X2] = ode45(@dinlin,[0 T],x0);
    end
 
    % grafico estados
    figure(1); clf;
    if ictr == 1
        subplot(3,2,1); plot(t/60,X(:,1),t2/60,X2(:,1)); xlabel('t (min)'); ylabel('V (m/s)');
        legend('sem PA','com PA');
        subplot(3,2,2); plot(t/60,X(:,2)*180/pi,t2/60,X2(:,2)*180/pi); xlabel('t (min)'); ylabel('gama (grau)');
        subplot(3,2,3); plot(t/60,X(:,3)*180/pi,t2/60,X2(:,3)*180/pi); xlabel('t (min)'); ylabel('alfa (grau)');
        subplot(3,2,4); plot(t/60,X(:,4)*180/pi,t2/60,X2(:,4)*180/pi); xlabel('t (min)'); ylabel('q (grau/s)');
        subplot(3,2,5); plot(t/60,X(:,5)/1000,t2/60,X2(:,5)/1000); xlabel('t (min)'); ylabel('H (km)');
        subplot(3,2,6); plot(X(:,6)/1000,X(:,5)/1000,X2(:,6)/1000,X2(:,5)/1000); xlabel('x (km)'); ylabel('H (km)');
    else
        subplot(3,2,1); plot(t/60,X(:,1),t2/60,X2(:,1)+Ve); xlabel('t (min)'); ylabel('V (m/s)');
        legend('nao linear','linear');
        subplot(3,2,2); plot(t/60,X(:,2)*180/pi,t2/60,X2(:,2)*180/pi); xlabel('t (min)'); ylabel('gama (grau)');
        subplot(3,2,3); plot(t/60,X(:,3)*180/pi,t2/60,(X2(:,3)+alfae)*180/pi); xlabel('t (min)'); ylabel('alfa (grau)');
        subplot(3,2,4); plot(t/60,X(:,4)*180/pi,t2/60,X2(:,4)*180/pi); xlabel('t (min)'); ylabel('q (grau/s)');
        subplot(3,2,5); plot(t/60,X(:,5)/1000,t2/60,(X2(:,5)+He)/1000); xlabel('t (min)'); ylabel('H (km)');
        subplot(3,2,6); plot(X(:,6)/1000,X(:,5)/1000); xlabel('x (km)'); ylabel('H (km)');
    end
     
    figure(2); clf;
    if ictr == 1
        subplot(1,2,1); plot(t2/60,ones(length(t2),1)*Ue(1)*180/pi,t2/60,(X2(:,9)+Ue(1))*180/pi); xlabel('t (min)'); ylabel('dp (grau)');
        subplot(1,2,2); plot(t2/60,ones(length(t2),1)*Ue(2)*100,t2/60,(X2(:,10)+Ue(2))*100); xlabel('t (min)'); ylabel('pi (%)');
        legend('sem PA','com PA');
    else
        subplot(1,2,1); plot(t/60,ones(length(t),1)*Ue(1)*180/pi); xlabel('t (min)'); ylabel('dp (grau)');
        subplot(1,2,2); plot(t/60,ones(length(t),1)*Ue(2)*100); xlabel('t (min)'); ylabel('pi (%)');
    end
 
    isim = input('Deseja simular outra perturbacao (sim=1)? ');
    disp ' '
     
end
 
%..........................................................................
 
% In Matlab:
% function dXdt = dinam(t,X,U);
% In Octave:
function dXdt = dinam(t,X);
global g S m c Iyy Cla Cld Cl0 Clq Clap Cm0 Cma Cmd Cmap Cmq Cd0 k...
    nv alfaf zf He Ve K ictr Xe Ue C Kr KI...
    U % In Octave
 
Vref = 240;
 
V = X(1); gama = X(2); alfa = X(3); q = X(4); H = X(5); x0 = X(6);
 
ro = atmosfera(H);
 
% PA
if ictr == 1
    w1 = X(7); w2 = X(8); ddp = X(9); dpif = X(10);
    dX = X(1:6) - Xe'; dX(6)=[];
    u  = -Kr*dX + KI*[w1;w2];
 
    dp = U(1) + ddp; pif = U(2) + dpif;
 
    ddpdt  = 10* (u(1)-ddp);
    dpifdt = 0.2*(u(2)-dpif);
 
    % batente
    if  dp  > 30*pi/180 & ddpdt > 0,  ddpdt = 0;  end
    if  dp  < -30*pi/180 & ddpdt < 0, ddpdt = 0;  end
    if  pif > 1 && dpifdt > 0,         dpifdt = 0; end
    if  pif < 0 && dpifdt > 0,         dpifdt = 0; end
    if abs(ddpdt) > 60*pi/180, ddpdt = sign(ddpdt)*60*pi/180; end
     
    dXdt(7,1)  = Ve - V;
    dXdt(8,1)  = He - H;
    dXdt(9,1)  = ddpdt;
    dXdt(10,1) = dpifdt;
else
    dp = U(1); pif = U(2);
end
 
 
Cl = Cl0 + Cla * alfa + Cld * dp + Clq*q*c/Vref;
Cd = Cd0 + k * Cl^2;
 
L = 1/2 * ro * V^2 * S * Cl;
D = 1/2 * ro * V^2 * S * Cd;
 
[Vwx Vwh dVwxdx0 dVwhdx0] = Vw(x0);
 
dXdt(6,1) = V*cos(gama) + Vwx;
 
dXdt(1,1) = (pif*Fmax(H,V)*cos(alfa + alfaf) - D - m*g*sin(gama)) / m ...
    - dVwxdx0*dXdt(6)*cos(gama) - dVwhdx0*dXdt(6)*sin(gama);
dXdt(2,1) = (pif*Fmax(H,V)*sin(alfa + alfaf) + L - m*g*cos(gama))/m/V ...
    + (dVwxdx0*dXdt(6)*sin(gama) - dVwhdx0*dXdt(6)*cos(gama))/V;
dXdt(3,1) = q - dXdt(2);
 
Cm = Cm0 + Cma*alfa + Cmd*dp + Cmap*dXdt(3)*c/Vref + Cmq*q*c/Vref;
M = 1/2 * ro * V^2 * S * c * Cm;
Mf = pif * Fmax(H,V) * zf;
 
dXdt(4,1) = (M + Mf)/Iyy;
dXdt(5,1) = V*sin(gama) + Vwh;
 
%..........................................................................
 
function f = equil(xe,ro)
global He Ve zf alfaf m g S c Cl0 Cla Cld Cd0 k Cm0 Cma Cmd
 
Fe = xe(1);
alfae = xe(2);
dpe = xe(3);
 
ro = atmosfera(He);
 
Cle = Cl0 + Cla * alfae + Cld * dpe;
Cde = Cd0 + k * Cle^2;
Cme = Cm0 + Cma*alfae + Cmd*dpe;
 
De = 1/2 * ro * Ve^2 * S * Cde;
Le = 1/2 * ro * Ve^2 * S * Cle;
Ma = 1/2 * ro * Ve^2 * S * c * Cme;
Mf = Fe * zf;
 
f1 = Fe*cos(alfae + alfaf) - De;
f2 = Fe*sin(alfae + alfaf) + Le - m*g;
f3 = Ma + Mf;
 
f = f1^2 + f2^2 + f3^2;
 
%..........................................................................
 
function dXdt = dinlin(t,X,U);
global A B
 
dXdt = A*X + B*U;
 
%..........................................................................
 
% vento
function [Vwx,Vwh,dVwxdx0,dVwhdx0] = Vw(x0)
global ax ah nx nh x0i x0f ivento
 
if x0 > x0i & x0 < x0f & ivento == 1
    Vwx     = ax                     * sin(2*pi*nx/(x0f-x0i)*(x0-x0f));
    Vwh     = ah                     * sin(2*pi*nh/(x0f-x0i)*(x0-x0f));
    dVwxdx0 = ax * 2*pi*nx/(x0f-x0i) * cos(2*pi*nx/(x0f-x0i)*(x0-x0f));
    dVwhdx0 = ah * 2*pi*nh/(x0f-x0i) * cos(2*pi*nh/(x0f-x0i)*(x0-x0f));
else
    Vwx     = 0;
    Vwh     = 0;
    dVwxdx0 = 0;
    dVwhdx0 = 0;
end
    
%..........................................................................
 
function F = eigfun(K,A,B,C)
F = sort(eig(A-B*K*C)); % Evaluate objectives