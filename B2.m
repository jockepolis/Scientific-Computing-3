%% Initial values, boundary conditions etc.
clear all; close all; clc;
%Choose number of grid points
N=input(' Choose a N = {41, 81, 161, 321, 641} : ');

h=1/(N-1);
T=1; %Final time
k=0.1*h;
t=0:k:T-k;
x=linspace(0,1,N);
u_0=zeros(1,N);
for i=1:N
    if abs(2*x(i)-0.3) <= 0.25
        u_0(i) = exp(-300*(2*x(i)-0.3)^2);
    elseif abs(2*x(i)-0.9) <= 0.2
        u_0(i) = 1;
    elseif abs(2*x(i)-1.6) <=0.2
        u_0(i) = sqrt(1-((2*x(i)-1.6)/0.2)^2);
    else
        u_0(i) = 0;
    end
end

%% Central difference discretization in space
u=zeros(N,round(T/k));
u(:,1)=u_0;
global A
A=diag(-ones(N-1,1),-1)+diag(ones(N-1,1),1);
A(1,end)=-1;
A(end,1)=1;
A=1/(2*h).*A;

%% Runge-Kutta 4 - discretization in time
F(round(T/k))=struct('cdata',[],'colormap',[]);
plot(x,u(:,1))
axis tight manual
ax=gca;
ax.NextPlot = 'replaceChildren';
anime=animatedline;
for i=1:round(T/k)-1
    K1=PDE(u(:,i),A);
    K2=PDE(u(:,i)+k/2*K1,A);
    K3=PDE(u(:,i)+k/2*K2,A);
    K4=PDE(u(:,i)+k*K3,A);
    u(:,i+1)=u(:,i)+k/6*(K1+2*K2+2*K3+K4);
    plot(x,u(:,i+1));
    title(['Solution at time= ',num2str(round(t(i+1),2))]);
    drawnow limitrate
    F(i+1)=getframe(gcf);
end
mesh(x,t,u');
xlabel('Space'); ylabel('Time'); title('Solution')

%% Error

%% Convergence

%% Functions
function dudt=PDE(u,A)
dudt=-1.*A*u;
end
