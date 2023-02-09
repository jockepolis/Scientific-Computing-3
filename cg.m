% Conjugate-gradient-solver function
function u=cg(A, b, TOL)
    u = zeros(size(b));
    r = A*u-b;
    d = r;
    rho=r'*r;
    k=0;
    while sqrt(rho)>TOL
        k = k+1;
        alpha= (rho)/(d'*A*d);
        u = u - alpha*d;
        r = r - alpha*A*d;
        rho_prev = rho; rho = r'*r;
        beta = rho/rho_prev;
        d = r + beta*d;
    end

end
