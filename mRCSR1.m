function [x, iter] = mRCSR1(f, x0, maxIter, tol)

%Parametros iniciales
    eta = 0.1;
    Delta = 1;
    maxDelta = 1.5;
    n = length(x0);
    r = 10^(-6); 
    iter = 0;
%Aproximacion del gradiente
    xk = x0;
    gk = apGrad(f, xk);
    Hk = speye(n);
    Bk = speye(n);
    
while norm(gk, inf) > tol && iter < maxIter
sk = -Hk*gk;
norms = norm(sk);
    if dot(sk, gk) < 0
        if  norms > Delta
            sk = ( Delta/norms )*sk;
        end
            
    else %Si -Hk*gk no es de decenso entonces tomamos el punto de Cauchy
        sk = pCauchy(Bk, gk, Delta);
    end
        
%Cociente de la reduccion en f y reduccion en m y ajuste del radio Delta
    redc = ( f(xk + sk) - f(xk) )/ ( dot(gk, sk) + 0.5 * sk'*Bk*sk);
    if redc < 0.1
        Delta = 0.5 * Delta;
        elseif redc > 0.75 &&  norm(sk) > 0.8 * Delta    
        Delta = min(maxDelta, 2*Delta);
    end
    
    if redc > eta
        xk = xk + sk;
        gplus = apGrad(f, xk);  
        gamma = (gplus - gk);
        gamma_Bs = gamma - Bk*sk;
        v = dot(gamma_Bs,sk);
        if abs(v) >= r * norm(sk) * norm(gamma_Bs) %Actualizamos para un r fijo.
            Bk = Bk + (1/v) * (gamma_Bs*gamma_Bs');
            w = sk - Hk*gamma;
            Hk = Hk + (1/dot(w, gamma)) * (w)*(w)';
        end
        gk = gplus;
        iter = iter + 1;
    end
        
end
    x = xk;   
end