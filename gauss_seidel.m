% Gauss-Seidel-solver function

% Set the starting conditions
function [itr_GS, x] = gauss_seidel(A, b, TOL)
x = zeros(size(b));
n=size(x,1);
norm=Inf;
itr = 0;
% The GS algorithm
while norm>TOL
    xold=x;
    for i=1:n
        y=0;
        for j=1:i-1
            y=y+A(i,j)*x(j);
        end
        for j=i+1:n
            y=y+A(i,j)*xold(j);
        end
        x(i)=(1/A(i,i))*(b(i)-y);
    end
    % Add another iteration unless norm < TOL, then break and print
    itr=itr+1;
    norm=abs(xold-x);
end
itr_GS = x;
fprintf('GS solves it in %f iterations \n',itr);
end
