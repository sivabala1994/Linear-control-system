
clc,clear
 m=1;
b=1;
r=1;
g=10;
J=1;
% CLTI System Matrice
 A = [0 1 0;0 -(b*r^2)/J (b*r)/J; -g (b*r)/m -b/m];
B = [0 0; 1 0 ; 0 1];
C=[1 0 0; 0 1 0; 0 0 1];

D = [0 0;0 0;0 0];

% Setting up a system
sysc = ss(A,B,C,D);

% setting up time
dt = .1; % choosing a time step size
tfin = 5;
tspan = 0:dt:tfin; % time span

% discretizing the system
sysd = c2d(sysc,dt); 

Ad = sysd.a;
Bd = sysd.b;
Cd = sysd.c;
Dd = sysd.d;

% Ad = expm(A*(tfin/dt));

% Checking for cotrollability and obervability
co = ctrb(Ad,Bd);
ob = obsv(Ad,Cd);

% Initial and final coditions
x0 = [pi/2 0 3]';
xf = [1 1 1]';

%% obtaining input to go from x0 to xf
co = [];
ob = [];

% costructing the cotrollability matrix
for k = 1:tfin/dt
    co = [co (A^((tfin/dt)-k))*B];
    ob = [ob;Cd*(A^(k-1))];
end

% Extracting the inputs
xbar = xf - (A^(tfin/dt))*x0;
[U1,S1,V1] = svd(co);
uco = V1*pinv(S1)*U1'*xbar;

% Separating the inputs

u1 = uco(1:2:length(uco));
u2 = uco(2:2:length(uco));
u1f=[u1,0];
u2f=[u2,0];
u=[u1,u2];
t=1:1:50
figure
plot(t,u','LineWidth',2),legend('input1','input2')
Simulation B*[u1(k);u2(k)];
xsc = x0;
for k = 1:length(u1)
    xsc(:,k+1) = Ad*xsc(:,k);
end

% Plotting
figure
plot(tspan,xsc,'LineWidth',2)

xlabel('time (s)')
ylabel('States')
legend('Attitude angle','angular velocity','lateral velocity')

initial condition

temp = zeros(size(C,1),1);
ybar = [];
for i = 1:tfin/dt
    for j = 1:i-1   % summation
        temp = temp + C*Ad^(i-1-j)*Bd* [u1(j);u2(j)];
    end
    cx = C*xsc(:,i);    % y(t) is c*x(t)
    ybar = [ybar; cx-temp];
    temp = zeros(size(C,1),1);
end

[U2,S2,V2] = svd(ob);
u_ini = V2*pinv(S2)*U2'*ybar;

closed loop
lambda=[-0.4 0.2 -0.6];
h=place(Ad,Bd,lambda);
xscc = x0;
for k = 1:length(u1)
    xscc(:,k+1) = (Ad-Bd*h)*xscc(:,k);
end

figure
plot(tspan,xscc,'LineWidth',2)

xlabel('time (s)')
ylabel('States')
legend('Attitude angle','angular velocity','lateral velocity')

observer
l2=[0.5;-0.01;-0.005]
L=place(Ad,C,l2);


yba=[];
for k = 1:length(u1)
    xscc(:,k+1) = (Ad-Bd*h)*xscc(:,k);
    yba(:,k+1)=C*xscc(:,k);

e(:,k+1)=xscc(:,k)-yba(:,k);
e_dot(:,k+1)=(Ad-C*L)*e(:,k);

end
figure
 plot(tspan,e_dot,'k-')
e=xscc-yba;
e_dot=(Ad-C*L)*e;

lambda=[-0.5 0.2 -0.6];
z=place(Ad,Bd,lambda);
xscn = [];
xscn(:,1)=[3.14,3.14,3];
xnl=[];
xnld=[];
xnl(:,1)=[3.14,3.14,3];
xnld(:,1)=[3.14,3.14,3];
for k = 1:length(u1)
    linear closed loop estimator/cpntroller
    xscn(:,k+1) = (Ad-Bd*z)*2*yba(:,k);
    
    
    cnl
    xnl(1,k+1)=xnl(2,k);
    xnl(2,k+1)=((b*r)/J)*yba(3,k)-((b*r^2)/J)*yba(2,k);
    xnl(3,k+1)=-(g)*sin(yba(1,k))- (b*r*yba(2,k))/m -(b*yba(3,k))/m;
    
   cnl with parameter variations
    xnld(1,k+1)=xnld(2,k);
    xnld(2,k+1)=(((b+1)*r)/J)*yba(3,k)-(((b+1)*r^2)/J)*yba(2,k);
    xnld(3,k+1)=-(g)*sin(yba(1,k))- ((b+1)*r*yba(2,k))/m -((b+1)*yba(3,k))/m;
    
end
figure
 plot(tspan,xnl(1,:),'k')


 
 
 
 hold on;
plot(tspan,xnld(1,:),'r')
legend('cnl','cnl with disturbance')

 figure
plot(tspan,e_dot,'LineWidth',2)
figure
 plot(tspan,xscc(3,:),'r')
 hold on;
 plot(tspan,e(3,:),'k')
legend('actual','estimated')
figure
 plot(tspan,xscc(2,:),'r')
 hold on;
 plot(tspan,e(2,:),'k')
legend('actual','estimated')
figure
 plot(tspan,xscc(1,:),'r')
 hold on;
 plot(tspan,e(1,:),'k')
legend('actual','estimated')

 figure
plot(tspan,xsc(1,:),'LineWidth',2)
hold on;
plot(tspan,xscc(1,:),'r','LineWidth',2)
xlabel('time')
ylabel('Angle')
legend('open loop','closed loop')
