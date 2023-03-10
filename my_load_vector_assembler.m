function B=my_load_vector_assembler(x)
%
% Returns the assembled load vector b.
% Input is a vector x of node coords.
%
fi = @(x) 2-tanh((x+1/2)/2*0.1);
N = length(x) - 1; B = zeros(N+1, 1); 
for i = 1:N 
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [fi(x(i)); fi(x(i+1))]*h/2;
end
B(end) = 7*1e+6; %to change the BC's