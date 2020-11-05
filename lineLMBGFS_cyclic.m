function [xk, iter] = lineLMBGFS_cyclic( f, x0, tol, maxiter, m )
% Purpose: approximate a local min of f using the linesearch algorithm
% and the (iBGFS) update formula (to avoid the solution of linear systems)
%
% Same parameters and results as lineDFP,
% with the exception that the (iBGFS) update formula is used.
%
    n = length(x0);
    xk = x0;
    g = apGrad(f,xk);
    H = speye(n);
    S = zeros(n,m);
    Gam = S;
    iter = maxiter;
    dk = -g;
    
    sA = 1;
    
    for k = 1:maxiter
        
        if norm(g, 'inf') <= tol
            iter = k-1;
            break
        end
        [alpha, gnew] = encAlpha( f, xk, dk, g );    %P2
        
        % memorizar
        s = alpha*dk;
        xk = xk + s;
        S(:,sA) = s;
        Gam(:,sA) = gnew-g;
        
        if k <m
            dk = -calcHg(S(:,1:k),Gam(:,1:k),gnew,sA);
        else
            dk = -calcHg(S,Gam,gnew,sA);
        end
        g = gnew;
        
        sA = mod(sA,m)+1;
        
    end
end

function [q] = calcHg(S,Gam, gnew,sA)
    q = gnew;
    m = length(S(1,:));
    irhos = 1 ./dot(S,Gam,1);
    alphas = zeros(m,1);
    indexes = mod((0:m-1)+ sA,m)+1;
    
    for i = fliplr(indexes)
        alphas(i) = dot(S(:,i), q)*irhos(i);
        q = q-alphas(i)*Gam(:,i);
    end
    delta = 1/(irhos(sA) * dot(Gam(:,sA),Gam(:,sA)));
    q = delta*q;
    
    for i = indexes
        beta = dot(Gam(:,i),q)*irhos(i);
        q = q + (alphas(i)-beta)*S(:,i);
        
    end
end