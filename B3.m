%% Initial values, boundary conditions etc.
clear all; close all; clc;
% Choose number of grid points
N=input(' Choose an N = {41, 81, 161, 321, 641} : ');

h=1/(N-1);
T=1; %Final time
k=0.1*h;
t=0:k:T-k;
x=linspace(0,1,N);
u_0=zeros(1,N);
% With pen and paper, Epsilon was determined to be Epsilon = h/2
% for the discretization to become a upwind method
epsilon = h/2;
%% Central difference discretization in space
u=zeros(N,round(T/k));
u(:,1)=u_0;
global A
A=diag((1/(2*h)+epsilon/h^2).*ones(N-1,1),-1)+diag((-2*epsilon/(h^2)).*ones(N,1));
A(1,end)=1/(2*h)+epsilon/h^2;

%% Forward Euler - discretization in time
for j=1:round(T/k)-1
    u(:,j+1) = u(:,j) + k*A*u(:,j);
end
