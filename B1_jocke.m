%% Initial values, boundary conditions etc.
clear all; close all; clc;
% Choose number of grid points
N=input(' Choose an N = {41, 81, 161, 321, 641} : ');

% Mesh in space
x=linspace(0,1,N);
h=1/(N-1);
% Mesh in time
T=1; %Final time
k=0.1*h;
t=0:k:T-k;
% BC
u_0=zeros(1,N);
for i=1:N
    if abs(2*x(i)-0.3) <= 0.25
        u_0(i)=exp(-300*(2*x(i)-0.3)^2);
    else
        u_0(i)=0;
    end
end
%% Central difference discretization in space
u=zeros(N,100);
u(:,1)=u_0;
global A
A=diag(-ones(N-1,1),-1)+diag(ones(N-1,1),1);
A(1,end)=-1;
A(end,1)=1;
A=1/(2*h).*A;

%% Runge-Kutta 4 - discretization in time
for i=1:T/k
    K1=PDE(u(:,i),A);
    K2=PDE(u(:,i)+k/2*K1,A);
    K3=PDE(u(:,i)+k/2*K2,A);
    K4=PDE(u(:,i)+k*K3,A);
    u(:,i+1)=u(:,i)+k/6*(K1+2*K2+2*K3+K4);
end
figure(1);
mesh(x,t,u');
xlabel('Space'); ylabel('Time')


%% Functions
function dudt=PDE(u,A)
dudt=-1.*A*u;
end
