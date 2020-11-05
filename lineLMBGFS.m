function [xk, iter] = lineLMBGFS( f, x0, tol, maxiter, m )
    n = length(x0);
    xk = x0;
    g = apGrad(f,xk);
    H = speye(n);
    S = zeros(n,m);
    Gam = S;
    iter = maxiter;
    dk = -g;
    for k = 1:maxiter
        if norm(g, 'inf') <= tol
            iter = k-1;
            break
        end
        [alpha, gnew] = encAlpha( f, xk, dk, g );    %P2
        
        % memorizar
        s = alpha*dk;
        xk = xk + s;
        S = [s, S(:,1:m-1)];
        Gam = [gnew-g, Gam(:,1:m-1)];
        
        if k <m
            dk = -calcHg(S(:,1:k),Gam(:,1:k),gnew);
        else
            dk = -calcHg(S,Gam,gnew);
        end
        g = gnew;
    end

end

function [q] = calcHg(S,Gam, gnew)
    q = gnew;
    m = length(S(1,:));
    irhos = 1./dot(S,Gam,1);
    alphas = zeros(m,1);
    for i = 1:m
        alphas(i) = dot(S(:,i), q)*irhos(i);
        q = q-alphas(i)*Gam(:,i);
    end
    
    delta = 1/(irhos(1)*dot(Gam(:,1),Gam(:,1)));
    q = delta*q;
    
    for i = m:-1:1
        beta = dot(Gam(:,i),q)*irhos(i);
        q = q + (alphas(i)-beta)*S(:,i);
        
    end
end