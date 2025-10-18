%% Linear Programming Control (LPC) Simulation for Rotary Inverted Pendulum
% This script implements the proposed LPC algorithm from the paper:
% "A Computationally Efficient Linear-Programming Control Framework for Constrained Input-Affine Nonlinear Systems"
% Author: Ali Safi, Ali Taghavian, Esmaeel Khanmirza, Fateme Namdarpour
% Repository: https://github.com/A-Safi/LPC
% 
% States: x = [alpha; alpha_dot; theta; theta_dot]
% Inputs : u = torque applied to the rotary arm (N·m)
%
% Dependencies on path:
%   - Controller.m
%   - IP.m           (signature: [fx,gx] = IP(x, params))
%   - RungeKutta.m   (signature: x_next = RungeKutta(@IP, params, x, t, u))
%
% The script demonstrates real-time variable-sampling LP-based control
% on the rotary inverted pendulum system.

clc;
clear;
% close all; %

% Add path to helper functions (Controller, IP, RungeKutta, etc.)
addpath("bin\")

% Output matrix (selects measurable variables: α and θ)
C = [0 1 0 0;
     0 0 0 1];

% Initial system states: [α; α̇; θ; θ̇]
x(:,1) = [0; 0; 0.1; 0];
time(1) = 0;   % initial time
U(1) = 0;      % initial control input
rd(1,:) = 0.1; % initial desired reference

% --- System and Control Parameters ---------------------------------------
params.m_pend = 0.0785;   % Pendulum mass (kg)
params.J_pend = 0.0007;   % Pendulum moment of inertia (kg·m²)
params.J_arm  = 0.0057;   % Arm moment of inertia (kg·m²)
params.l_pend = 0.25;     % Pendulum length (m)
params.r_arm  = 0.155;    % Arm length (m)
params.c_pend = 0.0006;   % Pendulum damping coefficient (N·s)
params.c_arm  = 0.035;    % Arm damping coefficient (N·s)
params.g      = 9.806;    % Gravitational acceleration (m/s²)

% --- Desired error-dynamics matrices ---
% These define how tracking errors decay under the control law
K1 = [-1, -10; 
       1.5, 14];
K2 = [-1.4, -1.4; 
       1.4, 1.4];

% --- Constraints on velocity, control input, and sampling period ---
dxmin = [-5; -5];   % Lower bounds on angular rates (rad/s)
dxmax = [ 5;  5];   % Upper bounds on angular rates (rad/s)
umin = -1.2;        % Minimum control torque (N·m)
umax =  1.2;        % Maximum control torque (N·m)
tmin = 0.01;        % Minimum sampling period (s)
tmax = 0.02;        % Maximum sampling period (s)

% --- Quick eigenvalue checks (stability verification) ---
eig(K1)               % Eigenvalues of K1 (desired dynamics)
eig(eye(2)+K1)        % Eigenvalues of I+K1
eig(inv(eye(2)+K2)*K1)% Eigenvalues of (I+K2)^(-1)*K1

% --- Reference signal (step command for α = π/4 after 100 samples) ---
ref = [0*ones(1,100), pi/4*ones(1,201)];

% --- Main Simulation Loop ------------------------------------------------
for k = 1:300
    % Compute nonlinear model terms f(x) and g(x)
    [fx, gx] = IP(x(:,k), params);

    % Desired post-update velocity (target for LPC)
    % Implements equation (8): ẋ⁺ = -K1·e - K2·ẋ
    dxplus = -K1 * [x(1,k) - ref(k); x(3,k)] ...
             -K2 * [x(2,k); x(4,k)];

    % Solve the linear programming problem
    % Controller() returns optimal decision vector p = [z; t; v]
    p = Controller(C*fx, C*gx, [x(2,k); x(4,k)], ...
                   dxmin, dxmax, umin, umax, ...
                   tmin, tmax, dxplus);
    
    % Recover true control input and sampling time
    v = p(end);        % product term (v = u*t)
    t = p(end-1);      % optimal sampling interval
    u = v / t;         % recover control input
    
    % Safety check for invalid control values
    if isnan(u)
        u = 0;
    end

    % Propagate system dynamics using variable-step Runge-Kutta integrator
    x(:,k+1) = RungeKutta(@IP, params, x(:,k), t, u);
    U(k+1) = u;           % store control input
    time_diff(k) = t;     % store sampling interval
    time(k+1) = time(k) + t; % update simulation time
end

%% Visualization
figure(1)

% α (pendulum angle)
subplot(3,2,1)
plot(time, x(1,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, ref, 'k-.', 'LineWidth', 1);
legend('Variable', 'Reference')
bnd = max(abs(x(1,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])
set(gca,'YTick', -pi/4:pi/4:pi/4)
set(gca,'YTickLabel', {'-π/4','0','π/4'})
ylabel('$\alpha$(rad)', 'interpreter','latex','FontSize',12,'FontName','Times')

% α̇ (pendulum angular velocity)
subplot(3,2,2)
plot(time, x(2,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, 3*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, -3*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
bnd = max(abs(x(2,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])
ylabel('$\dot\alpha$(rad/s)', 'interpreter','latex','FontSize',12,'FontName','Times')

% θ (arm angle)
subplot(3,2,3)
plot(time, x(3,:), 'b-', 'LineWidth', 1.2); hold on; grid on
ylabel('${\theta}$(rad)', 'interpreter','latex','FontSize',12,'FontName','Times')
bnd = max(abs(x(3,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])

% θ̇ (arm angular velocity)
subplot(3,2,4)
plot(time, x(4,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, 3*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, -3*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
bnd = max(abs(x(4,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])
ylabel('${\dot\theta}$(rad/s)', 'interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)', 'interpreter','latex','FontSize',12,'FontName','Times')

% Control input (torque)
subplot(3,2,5)
plot(time, U(1,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, 1.2*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, -1.2*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
ylabel('$T$(N.m)', 'interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)', 'interpreter','latex','FontSize',12,'FontName','Times')
bnd = max(abs(U(1,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])

% Variable sampling time evolution
figure
plot(time(1:end-1), time_diff, 'b*', 'MarkerSize', 5); axis tight; grid on; hold on
plot(time, tmin*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, tmax*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
ylabel('Sampling Time(sec)', 'interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)', 'interpreter','latex','FontSize',12,'FontName','Times')
ylim([0 0.025])