%% Linear Programming Control (LPC) Simulation for a Revolute–Prismatic (R–P) Robot
% This script implements the proposed LPC algorithm from the paper:
% "A Computationally Efficient Linear-Programming Control Framework for Constrained Input-Affine Nonlinear Systems"
% Author: Ali Safi, Ali Taghavian, Esmaeel Khanmirza, Fateme Namdarpour
% Repository: https://github.com/A-Safi/LPC
%
% States: x = [r; r_dot; theta; theta_dot]
% Inputs: u = [F; T]  (F = prismatic force (N), T = revolute torque (N·m))
%
% Dependencies on path:
%   - Controller.m
%   - RP.m           (signature: [fx,gx] = RP(x, params))
%   - RungeKutta.m   (signature: x_next = RungeKutta(@RP, params, x, t, u))

clc;
clear;
% close all;

addpath("bin\")  % path to Controller.m, RP.m, RungeKutta.m, etc.

% Output matrix (selects measurable variables: [r_dot; theta_dot])
C = [0 1 0 0;
     0 0 0 1];

% Initial Conditions -------------------------------
x(:,1)  = [0; 0; 0.5; 0];  % [r; r_dot; theta; theta_dot]
time(1) = 0;
U(:,1)  = [0; 0];


% ----------------------- R–P Model Parameters -----------------------------
params.a   = 0.5;
params.J1  = 0.0833;
params.J2  = 0.0017;
params.miu = 1;     
params.M   = 0.5;

% ----------------------- Controller Gains & Constraints -------------------
K1 = [2    0;
      0    1];
K2 = [0.05 0;
      0    0.01];

dxmin = [-0.3; -0.2];    % [r_dot_min; theta_dot_min]
dxmax = [ 0.3;  0.2];    % [r_dot_max; theta_dot_max]

umin  = [-0.5; -0.5];    % [F_min; T_min]
umax  = [ 0.5;  0.5];    % [F_max; T_max]

tmin  = 0.05;            % min sampling time (s)
tmax  = 0.1;             % max sampling time (s)

% Quick sanity checks (optional)
eig(K1)
eig(eye(2) + K2)
eig(inv(eye(2) + K2) * K1)
%%
% ----------------------- References ---------------------------------------
ref_r     = [zeros(1,100),  ones(1,200), -ones(1,500)];   % r (m)
ref_theta = [zeros(1,200), (pi/4)*ones(1,500)];           % theta (rad)

% ----------------------- Simulation ---------------------------------------
N = 550;

for k = 1:N
    % Nonlinear model (now with params)
    [fx, gx] = RP(x(:,k), params);

    % Desired post-update velocity: dxplus = -K1*e - K2*dx
    dxplus = -K1 * [x(1,k) - ref_r(k);
                    x(3,k) - ref_theta(k)] ...
             -K2 * [x(2,k);
                    x(4,k)];

    % Solve LP -> p = [z; t; v1; v2], where v = u * t
    p = Controller( C*fx, C*gx, [x(2,k); x(4,k)], ...
                    dxmin, dxmax, umin, umax, tmin, tmax, dxplus );

    % Recover sampling time and control input
    t = p(end-2);
    v = p(end-1:end);
    u = v ./ t;

    % Guard against infeasible/NaN outputs
    if any(isnan([u; t])) || t <= 0
        u = [0; 0];
        t = tmin;
        disp('LP infeasible — applying zero input for one step.');
    end

    % State propagation (pass params to the dynamics)
    x(:,k+1) = RungeKutta(@RP, params, x(:,k), t, u);

    % Logging
    U(:,k+1)     = u;
    time_diff(k) = t;
    time(k+1)    = time(k) + t;
end

%% ----------------------- Visualization -----------------------------------

figure(1)

% r (prismatic displacement)
subplot(2,2,1)
plot(time, x(1,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, ref_r(1:length(time)), 'k-.', 'LineWidth', 1);
ylabel('Displacement (m)', 'interpreter','latex','FontSize',12,'FontName','Times')
ub = max(x(1,:)); db = min(x(1,:) - 0.05);
ylim([1.1*db 1.1*ub]); xlim([0 time(end)])
legend('Variable', 'Reference')

% r_dot (prismatic velocity)
subplot(2,2,2)
plot(time, x(2,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, dxmax(1)*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, dxmin(1)*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
ylabel('Velocity (m/s)', 'interpreter','latex','FontSize',12,'FontName','Times')
ub = max(x(2,:)); db = min(x(2,:) - 0.05);
ylim([1.1*db 1.1*ub]); xlim([0 time(end)])

% theta (revolute angle)
subplot(2,2,3)
plot(time, x(3,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, ref_theta(1:length(time)), 'k-.', 'LineWidth', 1);
ylabel('Angle (rad)', 'interpreter','latex','FontSize',12,'FontName','Times')
xlabel('Time (sec)', 'interpreter','latex','FontSize',12,'FontName','Times')
ub = max(x(3,:)); db = min(x(3,:) - 0.05);
ylim([1.1*db 1.1*ub]); xlim([0 time(end)])
set(gca, 'YTick', -pi/2:pi/4:pi/2)
set(gca, 'YTickLabel', {'-π/2','-π/4','0','π/4','π/2'})

% theta_dot (revolute angular velocity)
subplot(2,2,4)
plot(time, x(4,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, dxmax(2)*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, dxmin(2)*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
ylabel('Angular Velocity (rad/s)', 'interpreter','latex','FontSize',12,'FontName','Times')
xlabel('Time (sec)', 'interpreter','latex','FontSize',12,'FontName','Times')
ub = max(x(4,:)); db = min(x(4,:) - 0.05);
ylim([1.1*db 1.1*ub]); xlim([0 time(end)])

figure(2)

% Input 1: prismatic force F
subplot(3,1,1)
plot(time, U(1,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, umax(1)*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, umin(1)*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
ylabel('$F$ (N)', 'interpreter','latex','FontSize',12,'FontName','Times')
ub = max(U(1,:)); db = min(U(1,:) - 0.1);
ylim([1.1*db 1.1*ub]); xlim([0 time(end)])

% Input 2: revolute torque T
subplot(3,1,2)
plot(time, U(2,:), 'b-', 'LineWidth', 1.2); hold on; grid on
plot(time, umax(2)*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, umin(2)*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
ylabel('$T$ (N$\cdot$m)', 'interpreter','latex','FontSize',12,'FontName','Times')
ub = max(U(2,:)); db = min(U(2,:) - 0.1);
ylim([1.1*db 1.1*ub]); xlim([0 time(end)])

% Variable sampling time
subplot(3,1,3)
plot(time(1:end-1), time_diff, 'b*', 'MarkerSize', 5); hold on; grid on; axis tight
plot(time, tmin*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
plot(time, tmax*ones(1,length(time)), 'r--', 'LineWidth', 1.2);
ylabel('Sampling Time (sec)', 'interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$ (sec)', 'interpreter','latex','FontSize',12,'FontName','Times')
ylim([0 0.15])
