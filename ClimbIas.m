function ClimbIas
  
  clear all
   
  global g = 9.80665;    % Gravity (m/s2)
  global m = 120000;     % Mass of the aircraft (kg)
  global S = 260;        % Area of reference (m2)
  global c = 6.61;       % Mean aerodynamic chord (m)
  global CLa =  4.982;   % Lift coefficient slope with angle of attack  
  global CLdp = .435;    % Lift coefficient slope with angle of horizontal stabilizer
  global CL0 = 0;        % Lift coefficient for zero angle of attack
  global Cm0 = -.025;    % Moment coefficient for zero angle of attack  
  global Cma = -1.246;   % Moment coefficient slope with angle of attack
  global Cmdp = -1.46;   % Moment coefficient slope with angle of horizontal stabilizer
  global alphaf = 1*pi/180; 
  global zf = 1.5;
  global Cnvi = 1.8e-5; 
  global mv = .62; 
  global mr = .122; 
  global Mi = .6; 
  global ri = .7768; 
  global roi = atmosfera(0);
  
  clc
  % Calculo do segmento de subida
  % Subida com Vias = cte e F = Fmax
  global Vias = 70;
  Hf   = 10000;
   
  % Simulaçao da subida
  X0    = [Vias 7*pi/180 0 m 0];
  [H,X] = ode45(@dinam,[0 Hf],X0);
   
  figure(1)
  subplot(2,2,1),plot(X(:,3)/60,X(:,1)); xlabel('t (min)'); ylabel('V (m/s)');
  subplot(2,2,2),plot(X(:,3)/60,X(:,2)*180/pi); xlabel('t (min)'); ylabel('gama (graus)');
  subplot(2,2,3),plot(X(:,5)/1000,H/1000); xlabel('X0 (km)'); ylabel('H (km)');
  subplot(2,2,4),plot(X(:,3)/60,X(:,4)); xlabel('t (min)'); ylabel('m (kg)');
endfunction
 
function dXdH = dinam(H,X);
  global zf alphaf g S c CL0 CLa CLdp ...
      Cm0 Cma Cmdp ...
      Cnvi mv mr Mi ri D...
      H2 V2 gama2 massa2 alpha
   
  V = X(1); gama = X(2); t = X(3); massa = X(4); x0 = X(5);
   
  [ro,drodh,Vsom] = atmosfera(H);
  M = V/Vsom;
   
  % alpha de subida
  H2 = H; V2 = V; gama2 = gama; massa2 = massa;

  % In Matlab:
  %alpha  = fzero(@alphafactor,7*pi/180);
  % In Octave:
  alpha  = sqp(7*pi/180,@alphafactor,[],[],-5*pi/180, 15*pi/180);
   
  F = Fmax(H,V);
   
  % deflexao do profundor - equil de momentos
  Mf  = F * cos(alphaf) * zf;
  Cmf = -Mf / (.5*ro*V^2*S*c);
  dp  = (Cmf-Cm0-Cma*alpha)/Cmdp;
   
  if round(alpha*10000)/10000 == 0.0882
      1
  end
  if isnan(alpha)
      1
  end
   
  CL = CL0 + CLa*alpha + CLdp*dp;
  CD = polar(CL,M);
  L  = .5*ro*V^2*S*CL;
  D  = .5*ro*V^2*S*CD;
   
  % consumo de combustivel
  Cnv = Cnvi * (M/Mi)^mv * (ro/ri)^mr;
  mcp = Cnv * F;
   
  dVdH    = (F*cos(alpha+alphaf) - D - massa*g*sin(gama))/massa   /(V*sin(gama));
  dgamadH = (F*sin(alpha+alphaf) + L - massa*g*cos(gama))/massa/V /(V*sin(gama));
  dtdH    = 1                                                     /(V*sin(gama));
  dmdH    = -mcp                                                  /(V*sin(gama));
  dx0dH   = V*cos(gama)                                           /(V*sin(gama));
   
  dXdH = [ dVdH ; dgamadH ; dtdH ; dmdH ; dx0dH ];
endfunction
 
function f = alphafactor(alpha)
  global roi S c Vias g alphaf zf H2 V2 gama2 massa2...
      CL0 CLa CLdp Cm0 Cma Cmdp
   
  [ro,drodh,Vsom] = atmosfera(H2);
  M = V2/Vsom;
   
  F = Fmax(H2,V2);
  Mf  = F * cos(alpha) * zf;
  Cmf = -Mf / (.5*ro*V2^2*S*c);
  dp  = (Cmf-Cm0-Cma*alpha)/Cmdp;
   
  CL = CL0 + CLa*alpha + CLdp*dp;
  CD = polar(CL,M);
  D  = 1/2*ro*V2^2*S*CD;
   
  %V = Vias * sqrt(roi/ro);
  Vp = -1/2 * Vias*sqrt(roi/ro^3) * drodh * V2*sin(gama2);
   
  f  = (F*cos(alpha + alphaf) - D - massa2*g*sin(gama2))/massa2 - Vp;
endfunction