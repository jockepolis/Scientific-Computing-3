% Temp.m
clear all
clf

% Model parameters
heatTransferCoefficient=1.63;
airTemperature=20;

% Domain discretization (Experiment these values!)
N=25;                   % the number of grid points
domainLength = 10;
dx=domainLength/(N+1); % calculate the grid size dx
x=[dx:dx:10-dx];       % grid points in space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we solve the system
%      LHS*u = rhs
% for the temperature solution u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic; % start the stopwatch timer

% Multiply both sides with dx^2
% Compute the right-hand-side vector rhs = -dx^2*h'*T_a
rhs=-heatTransferCoefficient*dx^2*airTemperature*ones(N,1);
% Impose boundary conditions
rhs(1)=rhs(1)-40;      % T = 40 at the left boundary
rhs(N)=rhs(N)-200;     % T = 200 at the right boundary

% Compute the matrix A approximating the second derivative multiplied
% by dx^2
A=zeros(N,N);
A(1,1)=-2;
A(1,2)=1;
for i=2:N-1
  A(i,i-1)=1;
  A(i,i)=-2;
  A(i,i+1)=1;
end
A(N,N-1)=1;
A(N,N)=-2;

% Compute the contribution of the reaction term -h'*T to the LHS matrix
R=zeros(N,N);
for i=1:N
  R(i,i)=-heatTransferCoefficient*dx^2;
end

% Assemble the left-hand-side matrix LHS
LHS=A+R;
% Solve for u
u=LHS\rhs;

t=toc; % stop the stopwatch timer

% Plot the solution
plot([0 x 10],[40 u' 200],'.-')
xlabel('Length')
ylabel('Temp')
disp(['Computational Time t = ' num2str(t) ' sec']);