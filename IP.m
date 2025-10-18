function [fx,gx]=IP(x,params)
%% Coefficients
mp = params.m_pend;
lp = params.l_pend;
jp = params.J_pend;
r  = params.r_arm;
ja = params.J_arm;  
c  = params.c_pend;
ca = params.c_arm;
gr = params.g;

a1 = (jp+mp*lp^2)*sin(x(3))^2 + mp*r^2 + ja;
a2 = mp*r*lp*cos(x(3));
a3 = 2*(jp+mp*lp^2)*x(4)*sin(x(3))*cos(x(3)) + ca;
a4 = mp*r*lp*x(4)*sin(x(3));
b1 = mp*r*lp*cos(x(3));
b2 = jp+mp*lp^2;
b3 = (jp+mp*lp^2)*x(2)*sin(x(3))*cos(x(3));
b4 = mp*gr*lp*sin(x(3));
c1 = 1 - (a2*b1/(a1*b2));

%% Outputs
fx = [x(2);
      1/(a1*c1)*(-(a2/b2)*(b4 + b3*x(2) - c*x(4)) - a3*x(2) + a4*x(4));
      x(4);
      1/(b2*c1)*( b4 - (b1/a1)*(-a3*x(2) + a4*x(4)) + b3*x(2) - c*x(4))];

gx = [0;
      1/(a1*c1);
      0;
     -1/(b2*c1)*(b1/a1)];
end
