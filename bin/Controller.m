function p = Controller(a,b,dxminus,dxa,dxb,umin,umax,tmin,tmax,dxhat)
%% Calculation
phi = dxhat-dxminus;

[n, m] = size(b) ;

M=[a,b];

H = [ M
     -M
     -umax, eye(m)
      umin,-eye(m)
      1, zeros(1,m)
     -1, zeros(1,m)];

g = [dxb-dxminus
     -dxa+dxminus
     zeros(m,1)
     zeros(m,1)
      tmax
     -tmin];

L = [-eye(n) -M
     -eye(n)  M];

q = [-phi
      phi];

H_p = [zeros(2*n+2*m+2,n), H];

Aineq = [L
         H_p];

bineq = [q
         g];

f = [ones(1,n), zeros(1,m+1)];

p = linprog(f,Aineq,bineq); % p_((n+1+m)*1) = [z_(n*1) delta_t_(1)  v_(m*1)]
