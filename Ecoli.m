% Ecoli.m
clear all
clf

% Problem parameters
D_D=0.28;
D_E=0.6;
sigma_1=20;
sigma_2=0.0063;
sigma_3=0.04;
sigma_4=0.8;
sigma_1p=0.028;
sigma_4p=0.027;

% Load initial state
load 'data/rho_DD.mat'
load 'data/rho_d.mat'
load 'data/rho_EE.mat'
load 'data/rho_e.mat'

% Space discretization
N=24;                 % number of grid points
dx=2/(N+1);           % space-step
x=[dx:dx:2-dx];       % grid points in space

% Time discretization
dt=0.0005;            % time step
Nt=200/dt;            % number of time steps

% Construct a matrix A approximating the second derivative
A=spdiags(repmat([1/dx^2 -2/dx^2 1/dx^2],N,1),[-1:1], N, N);
% To impose periodic boundary:
A(1,1) = -1/dx^2;
A(1,2) = 1/dx^2;
A(N,N-1) = 1/dx^2;
A(N,N) = -1/dx^2;

i=0;

tic; % start the stopwatch timer

for n=1:Nt
  if (mod(n,1000)==0)
    disp(['Timestep nr.',num2str(n),' of ',num2str(Nt)])
  end
  %
  % Approximate second order derivatives with D+D-
  %
  diff1 = D_D*A*rho_D;
  diff3 = D_E*A*rho_E;
  %
  % Evaluate right hand side terms
  %
  term1=sigma_1.*rho_D./(1.0+sigma_1p.*rho_e);
  term2=sigma_2.*rho_e.*rho_d;    
  term3=sigma_3.*rho_E.*rho_D;    
  term4=sigma_4.*rho_e./(1.0+sigma_4p.*rho_D);
  %
  % Time stepping with the Euler forward method
  %
  rho_D=rho_D+dt*(diff1-term1+term2);
  rho_d=rho_d+dt*(term1-term2);
  rho_E=rho_E+dt*(diff3-term3+term4);
  rho_e=rho_e+dt*(term3-term4);
  
  if (mod(n,10000)==0)
    i=i+1;
    rhosave_D(:,i)=rho_D+rho_d;
    rhosave_E(:,i)=rho_E+rho_e;
    t(i)=n*dt;
  end
end

comptime=toc; % stop the stopwatch timer

% Plot the figures
for k=1:4
  figure(k),clf
  set(gca,'FontSize',16)
end
  
[T X]=meshgrid(t,x);
figure(1)
surf(X,T,rhosave_D)
drawnow
view(2)
shading interp
H=colorbar; set(H,'FontSize',16)
xlabel('x')
ylabel('t')
title('\rho_D+\rho_d')

figure(2)
surf(X,T,rhosave_E)
drawnow
view(2)
shading interp
H=colorbar; set(H,'FontSize',16)
xlabel('x')
ylabel('t')
title('\rho_E+\rho_e')

figure(3)
plot(x,mean(rhosave_D'))
drawnow
xlabel('x')
title('Average of \rho_D+\rho_d')

figure(4)
plot(x,mean(rhosave_E'))
drawnow
xlabel('x')
title('Average of \rho_E+\rho_e')