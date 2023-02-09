function CG = testCG(A,b,TOL)
    x = b;
    r = b - A*x;
    if norm(r) < TOL
        return
    end
    y = -r;
    z = A*y;
    s = y'*z;
    t = (r'*y)/s;
    x = x + t*y;
  
    for k = 1:numel(b)
       r = r - t*z;
       if( norm(r) < TOL )
            return;
       end
       B = (r'*z)/s;
       y = -r + B*y;
       z = A*y;
       s = y'*z;
       t = (r'*y)/s;
       x = x + t*y;
    end
    CG = x;
 end