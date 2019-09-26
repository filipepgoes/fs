% Atmosphere
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