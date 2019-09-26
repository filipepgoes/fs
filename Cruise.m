function Cruise
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
  global alphaf = 1*pi/180; % Angle between the propulsion force and the longitudinal axis (rad)
  global zf = 1.5;          % Distance between the propulsion force and the longitudinal axis (m)
  global Cnvi = 1.8e-5;     % Coefficient of fuel consumption model
  global mv = .62;          % Mach number exponent for fuel consumption model
  global mr = .122;         % Density exponent for fuel consumption model
  global Mi = .6;           % Reference mach number for fuel consumption model
  global ri = .7768;        % Reference density for fuel consumption model
  global CLmax = 1;         % Maximum lift coefficient
  global qmax = 5.7e4;      % Maximum structural dynamic pressure 
  global sigma = .8;        % Coefficient of the optimization factor for the cruise conditions
  
  Vi=190;            % Initial guess for velocity (m/s2)
  Hi = 5000;                % Initial guess for altitude (m)
  clc

  % Initial conditions
  Xcruz = sqp([Hi Vi],@cruisefactor,[],@limits);
  Hi    = Xcruz(1);
  Vi    = Xcruz(2);

  % Cruise dynamics
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
  global He Ve alphaf m g S CL0 CLa CLdp...
      Cnvi mv mr Mi ri D Hcruz Vcruz

  V = X(1); gama = X(2); H = X(3); t = X(4); massa = X(5);

  [ro,drodh,Vsom] = atmosfera(H);
  M = V/Vsom;

  % Optimal conditions
  aux   = m;
  m     = massa;
  Xcruz = sqp([H V],@cruisefactor,[],@limits);%[5000 100]
  m     = aux;
  Hcruz = Xcruz(1);
  Vcruz = Xcruz(2);

  % Equilibrium
  aux    = m;
  m      = massa;
  xe0    = [m 0 0]; He = Hcruz; Ve = Vcruz;
  xe     = fminsearch(@equilibriumfactor,xe0);
  m      = aux;
  Fe     = xe(1);
  pife   = Fe/Fmax(He,Ve);
  alphae = xe(2);
  dpe    = xe(3);

  CL = CL0 + CLa*alphae + CLdp*dpe;
  CD = polar(CL,M);
  L  = 1/2 * ro * V^2 * S * CL;
  D  = 1/2 * ro * V^2 * S * CD;

  % Fuel consumption
  mcp = Fe * Cnvi * (M/Mi)^mv * (ro/ri)^mr;

  dVdx     = (Fe*cos(alphae+alphaf) - D - massa*g*sin(gama))/massa     /(V*cos(gama));
  dgamadx  = (Fe*sin(alphae+alphaf) + L - massa*g*cos(gama))/(massa*V) /(V*cos(gama));
  dHdx     = (V*sin(gama))                                             /(V*cos(gama));
  dtdx     = 1                                                         /(V*cos(gama));
  dmassadx = -mcp                                                      /(V*cos(gama));

  dfdx = [ dVdx ; dgamadx ; dHdx ; dtdx ; dmassadx ];
endfunction


function f = cruisefactor(X)
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

% Return array of constraints for drag (must be smaller than maximum thrust) and 
% velocity (greater than stall and lesser than maximum structural).
function c=limits(X)
  global Vs Vq D
  
  H = X(1); 
  V = X(2);
  
  c1 = Fmax(H,V) - D; % D < Fmax
  c2 = V - Vs;        % Vs < V
  c3 = Vq - V;        % V < Vq
  c = [c1 c2 c3];
endfunction


% Calculates a factor of difference between propulsive and resistive forces. 
% Minimizing it will give the control parameters for the desirable equilibrium 
% conditions.
function f = equilibriumfactor(xe)
  global He Ve zf alphaf m g S c CL0 CLa CLdp CD0 k Cm0 Cma Cmdp

  Fe = xe(1);     % Thrust 
  alphae = xe(2); % Angle of attack
  dpe = xe(3);    % Deflection of the horizonal stabilizer

  [ro,drodh,Vsom] = atmosfera(He);

  CLe = CL0 + CLa * alphae + CLdp * dpe;
  Me=Ve/Vsom;
  CDe = polar(CLe,Me);
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