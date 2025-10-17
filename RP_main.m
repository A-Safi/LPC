clc;
clear;
% close all;
addpath("bin\")
C=[0 1 0 0;0 0 0 1];

x(:,1) = [0; 0; 0.5; 0];
time(1) = 0;
U(:,1) = [0;0];
rd(1,:)=0.1;

tf2ss(1,conv([1 -1],[1 -2]));

K1 = [2 0;0 1];
K2 = [0.05, 0 ; 0, 0.01];
dxmin = [-0.3;-0.2];
dxmax = [0.3;0.2];
umin = [-0.5;-0.5];
umax = [0.5;0.5];
tmin = 0.05;
tmax = 0.1;

eig(K1)
eig(eye(2)+K2)
eig(inv(eye(2)+K2)*K1)

refp = [0*ones(1,100),1*ones(1,200),-1*ones(1,500)];
refr = [0*ones(1,200),pi/4*ones(1,500)];

for  k = 1:550
    [fx,gx] = RP(x(:,k));

    dxplus = -K1*[x(1,k)-refp(k);x(3,k)-refr(k)]-K2*[x(2,k);x(4,k)];
    
    p = Controller(C*fx,C*gx,[x(2,k);x(4,k)],dxmin,dxmax,umin,umax,tmin,tmax,dxplus);
    
    v = p(end-1:end);
    t = p(end-2);
    u = v./t;
    if isnan(u)
        u = [0;0];
        disp('not feas')
    end
    x(:,k+1) = RungeKutta(@RP,x(:,k)+0.0*[1;0;1;0].*rand(4,1)-0.0*[1;0;1;0],t,u);
    U(:,k+1) = u;
    time_diff(k) = t;
    time(k+1) = time(k) + t;
end

%%
figure(1)
subplot(2,2,1)
plot(time,x(1,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,refp(1:length(time)),'k-.','LineWidth',1);hold on;grid on
% title('Hi','FontSize',12,'FontName','Times')
ub=max(x(1,:));
db=min(x(1,:)-0.05);
ylim([1.1*db 1.1*ub])
xlim([0 time(end)])
ylabel('Displacement(m)','interpreter','latex','FontSize',12,'FontName','Times')


subplot(2,2,2)
plot(time,x(2,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,dxmax(1)*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,dxmin(1)*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('Velocity(m/s)','interpreter','latex','FontSize',12,'FontName','Times')
ub=max(x(2,:));
db=min(x(2,:)-0.05);
ylim([1.1*db 1.1*ub])
xlim([0 time(end)])

subplot(2,2,3)
plot(time,x(3,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,refr(1:length(time)),'k-.','LineWidth',1);hold on;grid on
ylabel('Angle(rad)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('Time(sec)','interpreter','latex','FontSize',12,'FontName','Times')
ub=max(x(3,:));
db=min(x(3,:)-0.05);
ylim([1.1*db 1.1*ub])
xlim([0 time(end)])
set(gca,'YTick',-pi/2:pi/4:pi/2) 
 set(gca,'YTickLabel',{'-π/2','-π/4','0','π/4','π/4'})
legend('Variable','Reference')

subplot(2,2,4)
plot(time,x(4,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,dxmax(2)*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,dxmin(2)*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('Angular Velocity (rad/s)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('Time(sec)','interpreter','latex','FontSize',12,'FontName','Times')
bnd=max(abs(x(4,:)));
ub=max(x(4,:));
db=min(x(4,:)-0.05);
ylim([1.1*db 1.1*ub])
xlim([0 time(end)])

figure(2)
subplot(3,1,1)
plot(time,U(1,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,umax(1)*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,umin(1)*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('$T$(N.m)','interpreter','latex','FontSize',12,'FontName','Times')
ub=max(U(1,:));
db=min(U(1,:)-0.1);
ylim([1.1*db 1.1*ub])
xlim([0 time(end)])

subplot(3,1,2)
plot(time,U(2,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,umax(2)*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,umin(2)*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('$T$(N.m)','interpreter','latex','FontSize',12,'FontName','Times')
ub=max(U(2,:));
db=min(U(2,:)-0.1);
ylim([1.1*db 1.1*ub])
xlim([0 time(end)])

subplot(3,1,3)
plot(time(1:end-1),time_diff,'b*','MarkerSize',5); axis tight;grid on;hold on
plot(time,tmin*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,tmax*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('Sampling Time(sec)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)','interpreter','latex','FontSize',12,'FontName','Times')
ylim([0 0.15])
