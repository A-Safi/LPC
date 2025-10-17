function [fx,gx]=RP(x)
%% Coefficients
a=0.5;
J1=0.0833;
J2=0.0017;
miu=1;
M=0.5;

%% Outputs
%d1 = 0*0.3*sin(2*((k-1)*dt))*exp(-0.2*((k-1)*dt));
%d2 = 0*0.4*cos(2*((k-1)*dt))*exp(-0.2*((k-1)*dt));
fx = [x(2);
    ((miu*x(1)+M*(x(1)+a))*(x(4)^2))/(miu+M)
    x(4)
    (-2*(miu*x(1)+M*(x(1)+a))*x(2)*x(4))/(J1+J2+miu*(x(1)^2)+M*((x(1)+a)^2))];
gx = [0 0;
    1/(miu+M) 0
    0 0
    0 1/(J1+J2+miu*(x(1)^2)+M*((x(1)+a)^2))];
end
