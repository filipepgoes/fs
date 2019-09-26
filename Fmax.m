% Maximum engine thrust
function y=Fmax(H,V)  
  
  nv = 0;        % Velocity exponent of the propulsive model
  nro=.75;       % Density exponent of the propulsive model
  roi = 1.225;   % Reference density of the propulsive model
  Fmaxi = 240000;% Maximum engine thrust (N)  
  Vtr=1;          % Base troposphere velocity (m/s2)
  Vst=1;          % Base stratosphere velocity
  
  ro=atmosfera(H);
  ros=atmosfera(11000);
  Fmaxs=Fmaxi*(ros/roi).^nro.*(Vst/Vtr).^nv;
  y1=(H<=11000).*(Fmaxi*(ro/roi).^nro.*(V/Vtr).^nv);
  y2=(H>11000).*(Fmaxs*(ro/ros).*(V/Vst).^nv);
  y=y1+y2;
endfunction