clc;
clear;
% close all;
addpath("bin\")
C=[0,1];

x(:,1) = [0; 0];
time(1) = 0;
U(1) = 0;
time1(1,:)=0;

%%
global m g l c
g=9.8;
m=0.5;
l=0.3;
c=0.1;

K = 5;
dxmin = -3;
dxmax = +3;
umin = -1.8;
umax = 1.8;
tmin = 0.01;
tmax = 0.01;

ref = [0*ones(1,50),pi*ones(1,251)];

for  k = 1:250
    [fx,gx] = SP(x(:,k));

    e = x(1,k) - ref(k);
    dxhat = -K * e;

    p = Controller(C*fx,C*gx,x(2,k),dxmin,dxmax,umin,umax,tmin,tmax,dxhat);
    
    v = p(end);
    t = p(end-1);%.01

    u =  v/t;%-(lambda * x(2,i) + K*(sign(s)) + fx(2))/gx(2);
    if isnan(u)
        u = 0;
    end
    
    x(:,k+1) = RungeKutta(@SP,x(:,k),t,u);
    U(k+1) = u;
    time_diff(k) = t;
    time(k+1) = time(k) + t;
end

%%
figure(1)
subplot(4,1,1)
plot(time,x(1,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,ref(1:length(time)),'k-.','LineWidth',1);hold on;grid on
% title('Hi','FontSize',12,'FontName','Times')
bnd=max(abs(x(1,:)));
ylim([-0.2 1.1*bnd])
xlim([0 time(end)])
legend('Variable','Refrence')
ylabel('$\theta$(rad)','interpreter','latex','FontSize',12,'FontName','Times')
set(gca,'YTick',-0:pi/4:pi) 
 set(gca,'YTickLabel',{'0','π/4','π/2','3π/4','π'})
subplot(4,1,2)
plot(time,x(2,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,3*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,-3*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('$\dot{\theta}$(rad/s)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)','interpreter','latex','FontSize',12,'FontName','Times')

bnd=max(abs(x(2,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])

subplot(4,1,3)
plot(time,U(1,:),'b-','LineWidth',1.2);hold on;grid on
plot(time,1.8*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,-1.8*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('$T$(N.m)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)','interpreter','latex','FontSize',12,'FontName','Times')

bnd=max(abs(U(1,:)));
ylim([-1.1*bnd 1.1*bnd])
xlim([0 time(end)])

subplot(4,1,4)
plot(time(1:end-1),time_diff,'b*','MarkerSize',5); axis tight;grid on;hold on
plot(time,tmin*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
plot(time,tmax*ones(1,length(time)),'r--','LineWidth',1.2);hold on;grid on
ylabel('Sampling Time(sec)','interpreter','latex','FontSize',12,'FontName','Times')
xlabel('$Time$(sec)','interpreter','latex','FontSize',12,'FontName','Times')
ylim([0 0.02])
