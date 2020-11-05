close; clear; clc;
N = [200, 1000];
M = [1,3,5,17,29];
tol = 1e-5;
maxiter = 2000;

f = @(x) dixmaan(x);

n = 200;
x_0 = 2*ones(n, 1);

tic
[~, iter] = lineBGFS(f, x_0, tol, maxiter);
time = toc;
fprintf('Metodo: LineBGFS,n: %d, iteraciones: %d, tiempo: %f\n',n,iter,time);

tic
[~, iter] = mRCSR1(f, x_0, maxiter, tol);
time = toc;
fprintf('Metodo: RCSR1,n: %d, iteraciones: %d, tiempo: %f\n',n,iter,time);
k = 1;
for n = N
    for m = M
        x_0 = 2*ones(n, 1);
        tic
        %[~, iter] = lineLMBGFS_cyclic(f, x_0, tol,maxiter, m);
         %Descomentar esta linea para usar el LMBGFS
        [~, iter] = lineLMBGFS(f, x_0, tol,maxiter, m);
        time = toc;
        fprintf('Metodo: LineBGFSLM Cyclic,n: %d,m: %d, iteraciones: %d, tiempo: %f\n',n,m,iter,time);
    end   
end
