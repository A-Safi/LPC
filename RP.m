function [fx,gx]=RP(x,params)
%% Coefficients
a   = params.a;
J1  = params.J1;
J2  = params.J2;
miu = params.miu;
M   = params.M;

%% Outputs
fx = [x(2);
    ((miu*x(1)+M*(x(1)+a))*(x(4)^2))/(miu+M)
    x(4)
    (-2*(miu*x(1)+M*(x(1)+a))*x(2)*x(4))/(J1+J2+miu*(x(1)^2)+M*((x(1)+a)^2))];
gx = [0 0;
    1/(miu+M) 0
    0 0
    0 1/(J1+J2+miu*(x(1)^2)+M*((x(1)+a)^2))];
end
