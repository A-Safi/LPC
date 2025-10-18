function p = Controller(a, b, dxminus, dxa, dxb, umin, umax, tmin, tmax, dxhat)
%% ------------------------------------------------------------------------
% Linear Programming Control (LPC) - Controller Function
%
% This function formulates and solves the linear optimization problem used
% in the LPC framework. It determines the optimal control input and sampling
% time that drive the system velocity (dx) toward the desired value (dxhat)
% while satisfying all system constraints.
%
% Inputs:
%   a, b      - System dynamics matrices (from f(x,xdot) and g(x,xdot))
%   dxminus   - Current state derivatives (velocities)
%   dxa, dxb  - Lower and upper bounds on dx
%   umin,
%   umax      - Control input bounds
%   tmin,
%   tmax      - Minimum and maximum allowable sampling times
%   dxhat     - Desired next-step derivative (target velocity)
%
% Output:
%   p         - Solution vector containing [z; delta_t; v]
%               where z is the absolute value slack, delta_t is the optimal
%               sampling time, and v = u * t (transformed variable)
%
% Author: Ali Safi
% ------------------------------------------------------------------------

%% Step 1. Compute the desired change in velocity
phi = dxhat - dxminus;

[n, m] = size(b);    % n: number of states, m: number of inputs

% Concatenate a and b for compact notation
M = [a, b];

%% Step 2. Define inequality constraints
% These inequalities enforce:
%   1) Velocity bounds (dxa ≤ dx ≤ dxb)
%   2) Input bounds (umin ≤ u ≤ umax)
%   3) Sampling time bounds (tmin ≤ t ≤ tmax)

H = [ M
     -M
     -umax, eye(m)
      umin,-eye(m)
      1, zeros(1,m)
     -1, zeros(1,m)];

g = [dxb - dxminus        % upper bound on dx
     -dxa + dxminus       % lower bound on dx
     zeros(m,1)           % upper bound on v (u*t) -> handled separately
     zeros(m,1)           % lower bound on v
     tmax                 % upper bound on sampling time
     -tmin];              % lower bound on sampling time

%% Step 3. Define L1-norm constraint representation
% Introduce auxiliary variable z to linearize |phi - (a*t + b*v)|
L = [-eye(n), -M          % z ≥  (phi - (a*t + b*v))
     -eye(n),  M];        % z ≥ -(phi - (a*t + b*v))

q = [-phi;                % right-hand side of inequalities
      phi];

%% Step 4. Assemble total constraint matrices for linprog
% z, delta_t, and v are optimization variables:
%   decision vector p = [z_(n×1); delta_t_(1×1); v_(m×1)]

% Add zero block for z terms in H constraints
H_p = [zeros(2*n + 2*m + 2, n), H];

% Combine all inequality constraints
Aineq = [L; H_p];
bineq = [q; g];

%% Step 5. Define cost function
% Objective: minimize sum(z)  →  minimize L1-norm of tracking error
f = [ones(1, n), zeros(1, m + 1)];

%% Step 6. Solve the Linear Program
% MATLAB’s built-in LP solver: minimize f*p  s.t. Aineq*p ≤ bineq
p = linprog(f, Aineq, bineq);

% The output vector p has structure:
%   p = [z_(n×1); delta_t_(1×1); v_(m×1)]
%
% From this, the actual control input is recovered as:
%   u = v / delta_t
%
% Notes:
% - This LP runs at every sampling step.
% - It enforces all constraints explicitly.
% - The solution ensures smooth, stable, and fast response.
% ------------------------------------------------------------------------
end
