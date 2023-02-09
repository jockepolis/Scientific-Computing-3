% PROBLEM A3
%%
clear all; clc; close all
format short
choice=input(' Choose an α = {1, 0.1, 0.001, 0.00001} : ');
if choice==1
    alpha=1;
elseif choice == 0.1
    alpha=0.1;
elseif choice == 0.001
    alpha=0.001;
elseif choice == 0.00001
    alpha=0.00001;
else
    choice=input('Try again, choose from {1, 0.1, 0.001, 0.00001} : ');
end
%%
n=10000;
A = matrix_machine(n,alpha);
b=randi(50,n,1);
TOL=10^(-5);
%%
format short
tic;
x1 = A\b ;
fprintf ('Backslash took % g sec \n ' , toc);
tic ;
x2 = jacobi (A, b, TOL);
fprintf ('Jacobi took % g sec \n ' , toc);
tic ;
x3 = gauss_seidel (A, b, TOL);
fprintf ('GS took % g sec \n ' , toc );
tic ;
x4 = cg (A, b, TOL);
fprintf ('CG took % g sec \n ' , toc );
tic ;
x5 = ownLUtest (A, b);
fprintf ('myownLU took % g sec \n ' , toc );
tic ;
[L, U] = lu(A);
y = L\b;
x6 = U\y;
fprintf ('Matlab LU took % g sec \n' , toc );

function A=matrix_machine(n,alpha)

alphas=(2+alpha)*ones(1,n);
A=diag(alphas)+diag(-1*ones(1,n-1),1)+diag(-1*ones(1,n-1),-1);

end