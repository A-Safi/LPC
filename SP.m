function [fx,gx]=SP(x,params)
%% Coefficients
m = params.m;
g = params.g;
l = params.l;
c = params.c;

%% Outputs
fx=[x(2);...
    -m*g*sin(x(1))/m/l-c*x(2)];
gx = [0;...
    1/m/(l^2);...
    ];
end

