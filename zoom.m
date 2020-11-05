function [alfa, gplus] = zoom(alfa_lo, alfa_hi, xk, dk,slope,f,fk)

c1 = 1e-4;
c2 = 0.99;
i = 1;
maxIter = 50;
F0 = f(xk + alfa_lo*dk);

while i <= maxIter
    alfa = (alfa_lo + alfa_hi)/2;
    Fi = f(xk + alfa*dk);
    gplus = apGrad(f, xk + alfa*dk);
    derivFi = dot(gplus, dk);
    L0 = fk + alfa*c1*slope;

    if ( Fi > L0 ) || ( Fi >= F0 )
        alfa_hi = alfa;
        i = i + 1;
        
    elseif (abs(derivFi) <= -c2*slope)
    break
    
    else
        
        if (derivFi*(alfa_hi - alfa_lo) >= 0)
            alfa_hi = alfa_lo;
        end
        
    alfa_lo = alfa;
    i = i + 1;
    F0 = Fi;
    end
end
end