function [X,er] = eitpMethod(A,X0,N,iterMax, tol)
% A: tensor de tamaño p x q x n
% alpha: parámetro escalar
% N: número de iteraciones RAPID internas (por defecto 5)
% tol: tolerancia para criterio de parada

% paper = Karmakar, B., & Behera, R. (2025). Efficient iterative methods for 
%         computing generalized inverse of tensors based on t-product: 
%         B. Karmakar, R. Behera. Computational and Applied Mathematics, 44(7), 380. 


% Dimensiones
[p, q, n] = size(A);

% Paso 2: Inicialización

% Paso 3: FFT a lo largo de la tercera dimensión
A_hat = fft(A, [], 3);
X_hat = fft(X0, [], 3);
I_hat = fft(eye(p), n, 3); % I: identidad extendida en 3D

Ainv=zeros(size(A,2),size(A,1),size(A,3));

% Paso 4: Iterar en cada frontal i
for i = 1:n
    X_hat_i=X_hat(:,:,i);
    A_hat_i=A_hat(:,:,i);
    I_hat_i=I_hat(:,:,i);
    for j=1:iterMax
        Pj = A_hat_i*X_hat_i;
        Uj = (1/4)*X_hat_i*(13*I_hat_i-Pj*(15*I_hat_i - Pj*(7*I_hat_i - Pj)));
        Vj = Uj + X_hat_i*(I_hat_i - A_hat_i*Uj) ;
        Yk=Uj;
        Wk=Vj;
        for k=1:N
            Zk=Wk+Yk*(I_hat_i-A_hat_i*Wk);
            Yk=Wk;
            Wk=Zk;
        end
        X_hat_i_New=Zk+ X_hat_i*(I_hat_i-A_hat_i*Zk);
        er=norm(X_hat_i_New-X_hat_i);
        if er<tol
            X_hat_i=X_hat_i_New;
            break
        end
        X_hat_i=X_hat_i_New;
    end
    Ainv(:,:,i)=X_hat_i;
end
X = ifft(Ainv, [], 3);
end