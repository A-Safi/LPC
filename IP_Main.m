clc;
clear;
% close all;
addpath("bin\")
C=[0 1 0 0;0 0 0 1];

x(:,1) = [0; 0; 0.1; 0];
time(1) = 0;
U(1) = 0;
rd(1,:)=0.1;
%%
global mp lp jp r ja c ca gr
mp = 0.0785;
jp = 0.0007;
ja = 0.0057;
lp = 0.25;
r  = 0.155;
c  = 0.0006;
ca = 0.035;
gr = 9.806;

K1 = [-1, -10 ; 1.5, 14];
K2 = [-1.4, -1.4 ; 1.4, 1.4];
dxmin = [-5;-5];
dxmax = [+5;+5];
umin = -1.2;
umax = 1.2;
tmin = 0.01;
tmax = 0.02;

eig(K1)
eig(eye(2)+K1)
eig(inv(eye(2)+K2)*K1)

ref = [0*ones(1,100),pi/4*ones(1,201)];

for  k = 1:300
    [fx,gx] = IP(x(:,k));

    dxplus = -K1*[x(1,k)-ref(k);x(3,k)]-K2*[x(2,k);x(4,k)];
    
    p = Controller(C*fx,C*gx,[x(2,k);x(4,k)],dxmin,dxmax,umin,umax,tmin,tmax,dxplus);
    
    v = p(end);
    t = p(end-1);
    u = v/t;
    if isnan(u)
        u = 0;
    end
    x(:,k+1) = RungeKutta(@IP,x(:,k),t,u);
    U(k+1) = u;
    time_diff(k) = t;
    time(k+1) = time(k) + t;
end

%%
figure(1)
subplot(3,2,1)
plot(time,x(1,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,ref,'k-.','LineWidth',1);hold on;grid on
legend('Variable','Reference')
% title('Hi','FontSize',12,'FontName','Times')
bnd=max(abs(x(1,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])
set(gca,'YTick',-pi/4:pi/4:pi/4) 
 set(gca,'YTickLabel',{'-π/4','0','π/4'})
ylabel('$\alpha$(rad)','interpreter','latex','FontSize',12,'FontName','Times')

subplot(3,2,2)
plot(time,x(2,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,3*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,-3*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
% title('Hi','FontSize',12,'FontName','Times')
bnd=max(abs(x(2,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])
ylabel('$\dot\alpha$(rad/s)','interpreter','latex','FontSize',12,'FontName','Times')

subplot(3,2,3)
plot(time,x(3,:),'b-','LineWidth',1.2);hold on;grid on
ylabel('${\theta}$(rad)','interpreter','latex','FontSize',12,'FontName','Times')
bnd=max(abs(x(3,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])

subplot(3,2,4)
plot(time,x(4,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,3*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,-3*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
% title('Hi','FontSize',12,'FontName','Times')
bnd=max(abs(x(4,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])
ylabel('${\dot\theta}$(rad/s)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)','interpreter','latex','FontSize',12,'FontName','Times')

subplot(3,2,5)
plot(time,U(1,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,1.2*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,-1.2*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('$T$(N.m)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)','interpreter','latex','FontSize',12,'FontName','Times')
bnd=max(abs(U(1,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])

figure
plot(time(1:end-1),time_diff,'b*','MarkerSize',5); axis tight;grid on;hold on
plot(time,tmin*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,tmax*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('Sampling Time(sec)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)','interpreter','latex','FontSize',12,'FontName','Times')
ylim([0 0.025])