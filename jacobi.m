% Jacobi-solver function

% Set the starting conditions
function itr_jacobi = jacobi(A, b, TOL)
x = zeros(size(b));
n=size(x,1);
norm=Inf;
itr = 0;
% The Jacobi algorithm
while norm>TOL
    xold=x;
    for i=1:n
        y=0;
        for j=1:n
            if j~=i
                y=y+A(i,j)*x(j);
            end
        end
        x(i)=(1/A(i,i))*(b(i)-y);
    end
    % Add another iteration unless norm < TOL, then break and print
    itr=itr+1;
    norm=abs(xold-x);
end
itr_jacobi = x;
fprintf('Jacobi solves it in %f iterations \n',itr);
end
