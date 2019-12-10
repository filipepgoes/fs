function arredondamento
clear all;
clc;

global Ap Bp Cp Hp x0...
    A B C H T ichute

% matrizes do sistema em malha aberta H=12m e V=70m/s
%       gama     alpha         q         H
Ap = [      0    0.4758   -0.0018   -0.0000
            0   -0.4758    1.0018    0.0000
            0   -0.6257   -0.2921   -0.0000
      70.0000         0         0         0];
Bp = [ 0.0404 ; -0.0404 ; -0.7709 ; 0];
Cp = [1 0 0 0];
Hp = [0 0 0 1];

% Sistema aumentado
A = [Ap Bp ; 0 0 0 0 -10];
B = [0;0;0;0;10];
C = [1 0 0 0  0];
H = [0 0 0 1  0];

% Kinicial - malha fechada
%        kp      ki      kd      kt      kr
K0   = [-0.2700 -0.0014 -0.1620 -1.8600  0.0360];
Kmin = -30*[1 1 1 1 1];
Kmax =  30*[1 1 1 1 1];
ichute = 1;
opt = optimset('MaxFunEvals',1000);
K0  = fmincon(@erro,K0,[],[],[],[],Kmin,Kmax,[],opt);

% simulink
BlocosArredondamento;
T = 20; % tempo de simulaçao
t = 0:.1:T;
r = [11*exp(-.3*t) ; -10.911*0.3048*exp(-.3*t)];

x0 = [-2.5*pi/180 0 0 16];
X0 = [x0 0 0];

ichute = 0;
K = fmincon(@erro,K0,[],[],[],[],Kmin,Kmax,[],opt);

kp = K(1); ki = K(2); kd = K(3);
kt = K(4); kr = K(5);

disp('ganhos encontrados');
disp(' kp ki kd kt kr');
disp(K);

% sistema com os ganhos definidos
ftr   = 1 / (1 + kd*H*B);
beta  = -ftr*(kp*H + kd*H*A + kt*C);
alpha = ftr*(kp + kr); 
gama  = ftr*ki;
delta = ftr*kd;
Ac    = [A+B*beta B*gama ; -H 0];
Bc    = [B*alpha B*delta ; 1 0]; 

disp ' '
disp('polos do sistema');
disp(eig(Ac));

% disp ' '
% disp(' Perturbaçoes no sistema:');
% dalpha = input('  Perturbaçao em alpha (grau): '); dalpha = dalpha*pi/180;
% dgama  = input('  Perturbaçao em gama (grau): '); dgama = dgama*pi/180;
% dq     = input('  Perturbaçao em q (grau/s): '); dq = dq*pi/180;
% dH     = input('  Perturbaçao em H (m): '); dH = dH*pi/180;
% 
% x0      = [-2.5*pi/180+dgama dalpha dq 16+dH 0 0];
% sist    = ss(Ac,Bc,[H 0],0);
% [y,t,x] = lsim(sist,r,t,x0);
% 
% figure(1);
% subplot(3,2,1); plot(t/60,x(:,1)*180/pi); xlabel('t (min)'); ylabel('gama (grau)');
% subplot(3,2,2); plot(t/60,x(:,2)*180/pi); xlabel('t (min)'); ylabel('alpha (grau)');    
% subplot(3,2,3); plot(t/60,x(:,3)*180/pi); xlabel('t (min)'); ylabel('q (grau/s)');
% subplot(3,2,4); plot(t/60,x(:,4));        xlabel('t (min)'); ylabel('H (m)');    
% subplot(3,2,5); plot(t/60,x(:,5));        xlabel('t (min)'); ylabel('dp (grau)');

function f = erro(K)
global Ap Bp Cp Hp x0...
    A B C H T ichute...
    %kp ki kd kt kr
 
kp = K(1); ki = K(2); kd = K(3);
kt = K(4); kr = K(5);

ftr  = 1 / (1 + kd*H*B);
beta = -ftr*(kp*H + kd*H*A + kt*C);
gama = ftr*ki;
Ac   = [A + B*beta B*gama ; -H 0];
rmax = max(real(eig(Ac)));

if ichute == 1, f = rmax;
else
    opt = simset('solver','ode45','SrcWorkspace','Current');
    [tout,xout,erro,dp] = sim('BlocosArredondamento',[0 T],opt);
    f = sum(abs(erro.*tout)) + 150*sum(abs(dp));
end
