% PROBLEM A2

w = 1;
N = 100;
TOL = 10e-5;
A = rand(N) + diag(w * ones(N, 1));
b = rand(N, 1);

tic;
x1 = A\b ;
fprintf ('Backslash took % g sec \n ' , toc);
tic ;
x2 = jacobi (A, b, TOL);
fprintf ('Jacobi took % g sec \n ' , toc);
tic ;
x3 = gauss_seidel (A, b, TOL);
fprintf ('GS took % g sec \n ' , toc );
% tic ;
% x4 = cg (A, b, TOL);
% fprintf ('CG took % g sec \n ' , toc );
tic ;
x5 = ownLUtest (A, b);
fprintf ('myownLU took % g sec \n ' , toc );
tic ;
[L, U] = lu(A);
y = L\b;
x6 = U\y;
fprintf ('Matlab LU took % g sec \n' , toc );