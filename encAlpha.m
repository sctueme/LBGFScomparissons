function [alfa, gplus] = encAlpha(f, xk, dk, gk)
    
c1 = 1e-4;
c2 = 0.99;
alfa0 = 0;
Maxalfa = 10;
i = 1;
maxIter = 50;
slope0 = dot(gk, dk);
fk = f(xk);
F0 = fk;

while i <= maxIter
    alfa = (alfa0+Maxalfa)/2; %Elegimos un alfa cualquiera en el intervalo
    Fi = f(xk + alfa*dk);
    gplus = apGrad(f, xk + alfa*dk);
    derivFi = dot(gplus, dk);
    L = fk + alfa*c1*slope0;
    
    %Elegimos alfa segun las condiciones de wolfe:
    if ( Fi > L ) || ( Fi >= F0 )
        [alfa, gplus] = zoom(alfa0, alfa, xk, dk, slope0, f, fk);
    break
        
    elseif (abs(derivFi) <= -c2*slope0)
    break
    
    elseif derivFi >= 0
    [alfa, gplus] = zoom(alfa, alfa0, xk, dk, slope0, f, fk);
    break
    
    else
    alfa0 = alfa;
    i = i + 1;
    F0 = Fi;
    end
end

    
end

