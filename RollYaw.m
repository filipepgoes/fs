function RollYaw
  clear all
  global CL0 CLa CLdp CLq CLap...
    CD0 k...
    Cm0 Cma Cmdp Cmq Cma Cmap...
    Cyb Cyda Cydd...
    Clb Clp Clr Clda Cldd...
    Cnb Cnp Cnr Cnda Cndd...
    zf alphaf nv nro...
    Ixx Iyy Izz Ixz...
    m g S c...
    He Ve betae psipte roe...
    ctr ti duracao...
    U % In Octave
 
  clc;
   
  CLa =  4.982; CLdp = .435; CL0 = 0; CLq = -.7; CLap = -.3;
  CD0 = 0.0175; k = 0.06;
  Cm0 = -.025; Cma = -1.246; Cmdp = -1.46; Cmap = -5; Cmq = -15;
  Cyb = -1.5; Cyda = .05; Cydd = .3;
  Clb = -1.3; Clp = -13; Clr = 2.9; Clda = -.33; Cldd = .25;
  Cnb =  1.75; Cnp = -1.5; Cnr = -7.5; Cnda = -.125; Cndd = -1;
  Fmax = 240000; nv = 0; nro =.75; alphaf = 1*pi/180; zf = 1.5;
  Ixx = 5.55e+6; Iyy = 9.72e+6; Izz = 14.51e+6; Ixz = -3.3e+4;
  m = 120000; g = 9.80665;    
  S = 260; c = 6.61;
   
   
  % entradas
  % He = input ('Altitude de voo (m): ');
  % Ve = input ('Velocidade de voo (m): ');
  % betae = input ('Angulo de derrapagem (\beta): ');
  % psipte = input ('Velocidade angular (\psi ponto): ');
   
      % equilibrio
          He     = 7000;
          Ve     = 200;
          betae  = 0*pi/180;
          psipte = 0*pi/180/60;
          roe    = atmosfera(He);
           
      % perturbaçao
          dV     = 0;         dH    = 200;
          dalpha = 0*pi/180;  dbeta = 2*pi/180;
          dtheta = 0*pi/180;  dfi   = 0*pi/180;
          dpsi   = 0*pi/180;  dp    = 0*pi/180;
          dq     = 0*pi/180;  dr    = 0*pi/180;
           
          dpif   = 0/100;     ddp   = 0*pi/180;
          da     = 0*pi/180;  dd    = 0*pi/180;
   
          ctr = 0; % 0=nada 1=redefine 2=pulso 3=double
          ctr = [ctr dpif ddp da dd];
          ti  = 1*60; duracao = 1;
           
          T = 10*60;
           
           
   
  % calculando o equilibrio
  %x = [alpha theta fi F dp dd da]
  x0 = [0 0 0 Fmax2(He,Ve) 0 0 0];
  xe = fminsearch(@equil,x0,optimset('MaxIter',50000,'MaxFunEvals',50000,'TolFun',1e-10,'TolX',1e-10));
   
  disp ' '
  disp('             NO EQUILIBRIO:');
  disp ' '
  disp('Tracao (N)'); Fe = xe(4); disp(Fe);
  disp('Posicao da manete (%)'); pife = Fe/Fmax2(He,Ve); disp(pife*100);
  disp('Angulo de ataque (graus)'); disp(xe(1)*180/pi); alphae = xe(1);
  disp('Angulo theta (graus)'); disp(xe(2)*180/pi); thetae = xe(2);
  disp('Angulo fi (graus)'); disp(xe(3)*180/pi); fie = xe(3);
  disp('Deflexao do profundor (graus)'); disp(xe(5)*180/pi); dpe = xe(5);
  disp('Deflexao do leme (graus)'); disp(xe(6)*180/pi); dde = xe(6);
  disp('Deflexao do aileron (graus)'); disp(xe(7)*180/pi); dae = xe(7);
   
  ue = Ve*cos(betae)*cos(alphae);
  ve = Ve*sin(betae);
  we = Ve*cos(betae)*sin(alphae);
   
  pe = -psipte*sin(thetae);
  qe =  psipte*cos(thetae)*sin(fie);
  re =  psipte*cos(thetae)*cos(fie);
   
  % linearizacao numerica
  Xe = [Ve alphae thetae qe He betae fie pe re 0 0 0];
  Ue = [pife ; dpe ; dae ; dde];
  eps1 = 1e-6;
  for i=1:9
      X = Xe; X(i) = Xe(i)+eps1; 
      % In Matlab:
      % f1 = dinam(0,X,Ue); 
      % In Octave:
      U=Ue; f1 = dinam(0,X); 
      f1(12)=[]; f1(11)=[]; f1(10)=[];
      X = Xe; X(i) = Xe(i)-eps1; 
      % In Matlab:
      % f2 = dinam(0,X,Ue); 
      % In Octave:
      U=Ue; f2 = dinam(0,X);
      f2(12)=[]; f2(11)=[]; f2(10)=[];
      A(:,i) = (f1 - f2) / (2*eps1);
  end
  for i=1:4
      U = Ue; U(i) = U(i)+eps1; 
      % In Matlab:
      % f1 = dinam(0,Xe,U); 
      % In Octave:
      f1 = dinam(0,Xe);
      f1(12)=[]; f1(11)=[]; f1(10)=[];
      U = Ue; U(i) = U(i)-eps1; 
      % In Matlab:
      % f2 = dinam(0,Xe,U); 
      % In Octave:
      f2 = dinam(0,Xe);
      f2(12)=[]; f2(11)=[]; f2(10)=[];
      B(:,i) = (f1 - f2) / (2*eps1);
  end
   
  disp ' '
  disp(' Movimento longitudinal')
  Along = A(1:5,1:5);
  Blong = B(1:2,1:2);
  damp(Along);
   
  disp ' '
  disp(' Movimento laterodirecional')
  Alatero = A(6:9,6:9);
  Blatero = B(6:9,3:4);
  damp(Alatero);
   
  % LQR
  [K,s,p] = lqr(A,B,eye(9),eye(4));
  disp('Ganhos encontrados:')
  disp(K)
  disp('Polos de malha fechada (LQR):')
  disp(eig(A-B*K))
   
  % simulaçao do equilibrio
  X0 = [Ve+dV alphae+dalpha thetae+dtheta qe+dq...
          He+dH betae+dbeta fie+dfi pe+dp...
          re+dr 0 0 0];
  U0 = [pife ; dpe ; dae ; dde];
  % In Matlab:
  % [t,X] = ode45(@dinam,[0 T],X0,[],U0);
  % In Octave:
  U = U0; [t,X] = ode45(@dinam,[0 T],X0);
  
  U = [0 0 0 0];
  for i = 1:length(t),
      uc = comandos(t(i));
      U(i,:) = uc + U0';
  end
   
  % X = [1  2      3      4  5  6     7  8 9 10 11 12];
  % X = [V  alfha  theta  q  H  beta  fi p r x0 y0 psi];
   
  % figura
  figure(1); clf;
  subplot(4,3,1);  plot(t/60,X(:,5)/1000);    xlabel('t (min)'); ylabel('H (km)');
  subplot(4,3,2);  plot(t/60,X(:,10)/1000);   xlabel('t (min)'); ylabel('x0 (km)');
  subplot(4,3,3);  plot(t/60,X(:,11)*180/pi); xlabel('t (min)'); ylabel('y0 (km)');
  subplot(4,3,4);  plot(t/60,X(:,1));         xlabel('t (min)'); ylabel('V (m/s)');
  subplot(4,3,5);  plot(t/60,X(:,2)*180/pi);  xlabel('t (min)'); ylabel('alpha (grau)');
  subplot(4,3,6);  plot(t/60,X(:,3)*180/pi);  xlabel('t (min)'); ylabel('theta (grau)');
  subplot(4,3,7);  plot(t/60,X(:,6)*180/pi);  xlabel('t (min)'); ylabel('beta (grau)');
  subplot(4,3,8);  plot(t/60,X(:,7)*180/pi);  xlabel('t (min)'); ylabel('fi (grau)');
  subplot(4,3,9);  plot(t/60,X(:,12)*180/pi); xlabel('t (min)'); ylabel('beta (grau)');
  subplot(4,3,10); plot(t/60,X(:,8)*180/pi);  xlabel('t (min)'); ylabel('p (grau/s)');
  subplot(4,3,11); plot(t/60,X(:,4)*180/pi);  xlabel('t (min)'); ylabel('q (grau/s)');
  subplot(4,3,12); plot(t/60,X(:,9)*180/pi);  xlabel('t (min)'); ylabel('r (grau/s)');
   
  figure(2); clf;
  subplot(3,2,1); plot3(X(:,10)/1000,X(:,11)/1000,X(:,5)/1000); xlabel('x0 (km)');   ylabel('y0 (km)'); zlabel('H (km)');
  subplot(3,2,2); plot(X(:,7)*180/pi,X(:,12)*180/pi);           xlabel('fi (grau)'); ylabel('beta (grau)');
  subplot(3,2,3); plot(t/60,U(:,1)*100);                        xlabel('t (min)');   ylabel('manete (%)');
  subplot(3,2,4); plot(t/60,U(:,2)*180/pi);                     xlabel('t (min)');   ylabel('profundor (grau)');
  subplot(3,2,5); plot(t/60,U(:,3)*180/pi);                     xlabel('t (min)');   ylabel('aileron (grau)');
  subplot(3,2,6); plot(t/60,U(:,4)*180/pi);                     xlabel('t (min)');   ylabel('leme (grau)');
endfunction 
% ______________________________________________________________
function uc = comandos(t);
  global ctr ti duracao
   
  if ctr(1) == 1
      if t > ti, uc = ctr(2:5);
      else uc = [0 0 0 0]; end
  elseif ctr(1) == 2
      if t > ti & t < ti+duracao, uc = ctr(2:5);
      else uc = [0 0 0 0]; end
  elseif ctr(1) == 3
      if t > ti & t < ti+duracao/2, uc = ctr(2:5);
      elseif t > ti+duracao/2 & t < ti+duracao, uc = -ctr(2:5);
      else uc = [0 0 0 0]; end
  elseif ctr(1) == 0, uc = [0 0 0 0]; end
endfunction
% ____________________________________________________________________
% In Matlab:
% function dXdt = dinam(t,X,U);
% In Octave:
function dXdt = dinam(t,X);
  
  global CL0 CLa CLdp CLq CLap...
      CD0 k...
      Cm0 Cma Cmdp Cmq Cma Cmap...
      Cyb Cyda Cydd...
      Clb Clp Clr Clda Cldd...
      Cnb Cnp Cnr Cnda Cndd...
      zf alphaf nv nro...
      Ixx Iyy Izz Ixz...
      m g S c...
      He Ve betae psipte roe...
      ctr ti duracao...
      U % In Octave 
  
  dXdt = zeros(9,1);
  Vref = 240;
   
  V  = X(1); alpha = X(2);  theta = X(3);  q   = X(4);
  H  = X(5); beta  = X(6);  fi    = X(7);  p   = X(8);
  r  = X(9); x0    = X(10); y0    = X(11); psi = X(12);
   
  uc = comandos(t);
   
  pif = U(1) + uc(1); dp = U(2) + uc(2);
  da  = U(3) + uc(3); dd = U(4) + uc(4);
   
  % if t ~= tanterior
  %     controles(contador,:) = [t pif dp da dd];
  %     contador = contador+1;
  %     tanterior = t;
  % end
   
  ro = atmosfera(H);
   
  u = V*cos(beta)*cos(alpha);
  v = V*sin(beta);
  w = V*cos(beta)*sin(alpha);
   
  CL = CL0 + CLa*alpha + CLdp*dp * CLq*q*c/V;
  CD = CD0 + k*CL^2;
  Cy = Cyb*beta + Cyda*da + Cydd*dd;
   
  CX = - CD*cos(alpha)*cos(beta) - Cy*cos(alpha)*sin(beta) + CL*sin(alpha);
  CY =                             Cy                                         ;
  CZ = - CD*sin(alpha)*cos(beta) - Cy*sin(alpha)*sin(beta) - CL*cos(alpha);
   
  X = 1/2 * ro * V^2 * S * CX;
  Y = 1/2 * ro * V^2 * S * CY;
  Z = 1/2 * ro * V^2 * S * CZ;
   
  Cl = Clb*beta + Clp*p*c/V + Clr*r*c/V + Clda*da + Cldd*dd;
  Cm = Cm0 + Cma*alpha + Cmdp*dp + Cmq*q*c/V;
  Cn = Cnb*beta + Cnp*p*c/V + Cnr*r*c/V + Cnda*da + Cndd*dd;
   
  L = 1/2 * ro * V^2 * S * c * Cl;
  M = 1/2 * ro * V^2 * S * c * Cm;
  N = 1/2 * ro * V^2 * S * c * Cn;
   
  F = pif * Fmax2(H,V);
  Mf = F * cos(alphaf) * zf;
       
  up = (- m*(w*q-v*r) - m*g*sin(theta)          + X + F*cos(alphaf))/m;
  vp = (- m*(u*r-w*p) + m*g*cos(theta)*sin(fi)  + Y                )/m;
  wp = (- m*(v*p-u*q) + m*g*cos(theta)*cos(fi)  + Z - F*sin(alphaf))/m;
  Vp = (u*up + v*vp + w*wp)/V;
   
  betap  = (V*vp-v*Vp) / (V^2*cos(beta));
  alphap = (u*wp - w*up) / (u^2 + w^2);
       
  pp = (-Izz*(L + Ixz*p*q) - Ixz*(N + (Ixx - Iyy)*p*q) + (Ixz^2 + Izz*(-Iyy + Izz))*q*r)/(Ixz^2 - Ixx*Izz);
  qp = (-(Ixx - Izz)*p*r - Ixz*(p^2 - r^2) + M + Mf)/Iyy;
  rp = -(Ixz^2*p*q + Ixx*(N + (Ixx - Iyy)*p*q) + Ixz*(L - (Ixx - Iyy + Izz)*q*r))/(Ixz^2 - Ixx*Izz);
   
  psip   =                  (q*sin(fi) + r*cos(fi))/cos(theta);
  thetap =                   q*cos(fi) - r*sin(fi);
  fip    = p + tan(theta) * (q*sin(fi) + r*cos(fi));
   
  Hp  = u*sin(theta)          - v*cos(theta)*sin(fi)                               - w*cos(theta)*cos(fi);
  x0p = u*cos(theta)*cos(psi) + v*(sin(fi)*sin(theta)*cos(psi) - cos(fi)*sin(psi)) + w*(cos(fi)*sin(theta)*cos(psi) + sin(fi)*sin(psi));
  y0p = u*cos(theta)*sin(psi) + v*(sin(fi)*sin(theta)*sin(psi) + cos(fi)*cos(psi)) + w*(cos(fi)*sin(theta)*sin(psi) - sin(fi)*cos(psi));
   
      dXdt(1)  = Vp;     dXdt(2)  = alphap;
      dXdt(3)  = thetap; dXdt(4)  = qp;
      dXdt(5)  = Hp;     dXdt(6)  = betap;
      dXdt(7)  = fip;    dXdt(8)  = pp;
      dXdt(9)  = rp;     dXdt(10) = x0p;
      dXdt(11) = y0p;    dXdt(12) = psip;
endfunction 
% ____________________________________________________________________
function f = equil(xe)
  global CL0 CLa CLdp CLq CLap...
      CD0 k...
      Cm0 Cma Cmdp Cmq Cma Cmap...
      Cyb Cyda Cydd...
      Clb Clp Clr Clda Cldd...
      Cnb Cnp Cnr Cnda Cndd...
      zf alphaf nv nro...
      Ixx Iyy Izz Ixz...
      m g S c...
      He Ve betae psipte roe
   
  alphae = xe(1);
  thetae = xe(2);
  fie    = xe(3);
  Fe     = xe(4);
  dpe    = xe(5);
  dde    = xe(6);
  dae    = xe(7);
         
  ue = Ve*cos(betae)*cos(alphae);
  ve = Ve*sin(betae);
  we = Ve*cos(betae)*sin(alphae);
   
  pe = - psipte*sin(thetae);
  qe =   psipte*cos(thetae)*sin(fie);
  re =   psipte*cos(thetae)*cos(fie);
   
  CLe = CL0 + CLa*alphae + CLdp*dpe * CLq*qe*c/Ve;
  CDe = CD0 + k*CLe^2;
  Cye = Cyb*betae + Cyda*dae + Cydd*dde;
   
  Cle = Clb*betae + Clp*pe*c/Ve + Clr*re*c/Ve + Clda*dae + Cldd*dde;
  Cme = Cm0 + Cma*alphae + Cmdp*dpe + Cmq*qe*c/Ve;
  Cne = Cnb*betae + Cnp*pe*c/Ve + Cnr*re*c/Ve + Cnda*dae + Cndd*dde;
   
  Cx = - CDe*cos(alphae)*cos(betae) - Cye*cos(alphae)*sin(betae) + CLe*sin(alphae);
  Cy =                                Cye                                         ;
  Cz = - CDe*sin(alphae)*cos(betae) - Cye*sin(alphae)*sin(betae) - CLe*cos(alphae);
   
  X = 1/2 * roe * Ve^2 * S * Cx;
  Y = 1/2 * roe * Ve^2 * S * Cy;
  Z = 1/2 * roe * Ve^2 * S * Cz;
  Mf = Fe*cos(alphaf)*zf;
   
  % força
  f1 = - m*(qe*we-ve*re) - m*g*sin(thetae)          + X + Fe*cos(alphaf);
  f2 = - m*(ue*re-pe*we) + m*g*cos(thetae)*sin(fie) + Y                 ;
  f3 = - m*(pe*ve-ue*qe) + m*g*cos(thetae)*cos(fie) + Z - Fe*sin(alphaf);
   
  % momento
  f4 = - (Izz-Iyy)*qe*re + Ixz*pe*qe       + 1/2*roe*Ve^2*S*c*Cle     ;
  f5 = - (Ixx-Izz)*pe*re - Ixz*(pe^2-re^2) + 1/2*roe*Ve^2*S*c*Cme + Mf;
  f6 = - (Iyy-Ixx)*pe*qe - Ixz*re*qe       + 1/2*roe*Ve^2*S*c*Cne     ;
   
  % altitude (H ponto)
  f7 = ue*sin(thetae) - ve*cos(thetae)*sin(fie) - we*cos(thetae)*cos(fie);
   
  f = f1^2 + f2^2 + f3^2 + f4^2 + f5^2 + f6^2 + f7^2;
endfunction 

