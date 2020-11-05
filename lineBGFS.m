function [xk, iter] = lineBGFS( f, x0, tol, maxiter )
    n = length(x0);
    xk = x0;
    g = apGrad(f,xk);
    H = speye(n);
    iter = maxiter;
    for k = 1:maxiter
        if norm(g, 'inf') <= tol
            iter = k-1;
            break
        end
        dk = -H*g;
        [alpha, gnew] = encAlpha( f, xk, dk, g );
        s = alpha*dk;
        xk = xk + s;
        gamma = gnew -g;
        irho = dot(gamma, s);
        Hgam = H*gamma/irho;
        H = H - (s*Hgam' + Hgam*s') + ((dot(gamma,Hgam) + 1)/irho*s)*s';
        g = gnew;
    end
end