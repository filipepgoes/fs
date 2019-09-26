% Drag polar
function CD=polar(CL,M)
  CD0 = 0.0175;   % Drag coefficient for zero lift coefficient angle of attack
  k = 0.06;       % Slope of drag coefficient with the square of lift coefficient
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