function [fx,gx]=SP(x)
%% Coefficients
global m  g l c

%% Outputs
fx=[x(2);...
    -m*g*sin(x(1))/m/l-c*x(2)];
gx = [0;...
    1/m/(l^2);...
    ];
end

