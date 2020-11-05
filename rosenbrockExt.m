function [resp] = rosenbrockExt(x,c)
    n = length(x);
    sum = 0;
    
    for i = 2:2:n
        xnext = x(i);
        xi = x(i-1);
        sum = sum + c*(xnext - xi^2)^2 + (1 - xi)^2;
    end
    resp = sum;
end