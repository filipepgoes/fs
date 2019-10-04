% ____________________________________________________________________
function y = Fmax2(H,V)  
  global nv nro
   
  if(H <= 11000)
      nro2 = nro;
  else
      nro2 = 1;
  end
   
  % ao nivel do mar
  roi = atmosfera(0);
  Vi = 120;
  fmaxi = 240000;
   
  % na tropopausa
  ros = atmosfera(11000);
  Vs = 120;
  fmaxs = fmaxi * (Vs/Vi).^nv * (ros/roi).^nro2;
   
  ro = atmosfera(H);
  [jV,jH] = size(H);
   
  if(H <= 11000)
      y = fmaxi * (V/Vi).^nv * (ro/roi).^nro2;
  else
      y = fmaxs * (V/Vs).^nv * (ro/ros);
  end
endfunction 