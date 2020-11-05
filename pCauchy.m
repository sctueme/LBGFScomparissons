function [pC] = pCauchy( B, g, delta )
% In : B ... (symmetric matrix) approximates the hessian of f in xk
%       g     ... (vector) gradient of  f  in  xk
%       delta ... trust region radius
% Out:  pC    ... The  Cauchy  point
    pk = g/norm(g);
    pkBpk = dot(pk,B*pk);
    if(pkBpk <= 0)
        t = delta;
    else
        t = min(delta,norm(g)^3/(delta*pkBpk));
    end
    pC = -0.99*t*pk;
end