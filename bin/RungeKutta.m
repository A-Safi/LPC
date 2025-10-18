function nstate = RungeKutta(fnc,params,state,h,u)

[fx,gx] = fnc(state,params);
k1 = fx+gx*u;

[fx,gx] = fnc(state+h/2*k1,params);
k2 = fx+gx*u;

[fx,gx] = fnc(state+h/2*k2,params);
k3 = fx+gx*u;

[fx,gx] = fnc(state+h*k3,params);
k4 = fx+gx*u;

nstate=state+h/6*(k1+2*k2+2*k3+k4);

end

