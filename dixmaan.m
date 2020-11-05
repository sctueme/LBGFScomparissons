function [resp] = dixmaan(x)
    % Using K Parameters
    % K2, K3 = 0
    alpha = 1;
    beta = 0.125;
    gamma = 0.125;
    delta = 0.125;
    k1 = 2;
    k4 = 2;
    n = length(x);
    m = floor(n/3);
    fs = 0;
    ss = 0;
    ts = 0;
    qs = 0;
    for i=1:n
        fs = fs + alpha*x(i)^2*(i/n)^k1;
    end
    for i=1:n-1
        ss = ss + beta*x(i)^2*(x(i+1)+x(i+1)^2)^2*1;
    end
    for i=1:2*m
        ts = ts + gamma*x(i)^2*x(i+m)^4*1;
    end
    for i=1:m
        qs = qs + delta*x(i)*x(i+2*m)*(i/n)^k4;
    end
    resp = 1 + fs +ss +ts +qs;
end