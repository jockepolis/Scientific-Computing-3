% Temperature.m
clear all
clf

% Model parameters
rho=2300;                      % Density
c=800;                         % Specific heat 
kappa=1.63;                    % Heat conduction 
T=3600;                        % Solution time

% Domain discretization (Experiment these values!)
N=151;                          % number of grid points
dt=0.2;                          % time step

nt=round(T/dt); % Number of steps
dx=0.2/(N+1);                  % space-step
x=[dx:dx:0.2-dx];              % grid points in space

t=0;                           % initiate time
u(1:N,1)=20*ones(N,1);         % initial solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize the equation with with central difference in space and
% forward Euler in time
% Read more: https://en.wikipedia.org/wiki/FTCS_scheme
% We now construct a matrix A such that the solution is updated as
% u = A * u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% setup the tridiagonal coefficient-matrix for the right-hand side %%%
c1=dt*kappa/(rho*c*dx^2);  %precompute common factor
A=spdiags(repmat([c1 1-2*c1 c1],N,1),[-1:1], N, N); %form the matrix

for n=1:nt
  u(:)=A*u(:);                  %compute right-hand side
  u(1)=u(1)+c1*20;              %set boundary condition at x=0
  u(N)=u(N)+c1*20-c1*50*t/3600; %set boundary condition at x=0.2 
  t=t+dt;
end

%plot results
clf;
subplot(2,1,1);
plot([0 x 0.2],[20 u' 20-50*t/3600],'bo:')
xlabel('x')
ylabel('Temperature')
title('Temperature at t=3600s')

%load precomputed reference solution
load data/U.mat;
load data/X.mat;
hold on;
plot([0 X' 0.2],[20 U' 20-50*t/3600],'r-')
legend('Numerical Solution','Reference Solution')


% Error
subplot(2,1,2)
y=spline(X,U,x)';
plot([0 x 0.2],[0 (y-u)' 0],'go:');axis([0 0.2 min([y-u;0]) max([y-u;0])]);
xlabel('x')
ylabel('Error')
title('Error in the computed solution at t=3600s')
