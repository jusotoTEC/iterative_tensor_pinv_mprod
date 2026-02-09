function Xk = cgtpMethod(A, X0, iterMax,tol)

% This function estimates the pseudoinverse of a tensor by computing
% the tensor A.
%
% References: B. Huang, Conjugate gradient-type method for the tensor linear
%             system via the T-product and its application in the calculation
%             of the Moore-Penrose inverse, Applied Mathematics and
%             Computation, vol. 472, p. 128627, 2024.
%
% Inputs: tensor A of size m x n x p.
%         Maximum number of iterations (MaxIter)
%         Tolerance (tol)
%
% Outputs: pseudoinverse Xk of size n x m x p
%
% Implemented by Samuel Valverde Sánchez and Juan Pablo Soto Quirós

    [m,n,p] = size(A);
    Y0 = zeros(m,n,p);
    AT = tCTranspose(A);
    R0_1 = A - tprod3(A,X0,A); R0_2 = -X0 + tprod3(AT,Y0,AT);
    P0_1 = tprod3(AT,R0_1,AT)+R0_2; P0_2 = -tprod3(A,R0_2,A);
    Q0_1 = P0_1; Q0_2 = P0_2;
    k = 0;
    while k < iterMax
        ak = (tNorm2(R0_1)^2 + tNorm2(R0_2)^2)/(tNorm2(Q0_1)^2 + tNorm2(Q0_2)^2);
        X = X0 + ak*Q0_1;
        Y = Y0 + ak*Q0_2;
        R1_1 = A - tprod3(A,X,A);
        R1_2 = -X + tprod3(AT,Y,AT);
        if tNorm2(R1_1)^2 + tNorm2(R1_2)^2 < tol
            break
        end
        P1_1 = tprod3(AT,R1_1,AT)+R1_2;
        P1_2 = -tprod3(A,R1_2,A);
        bk = (tNorm2(R1_1)^2 + tNorm2(R1_2)^2)/(tNorm2(R0_1)^2 + tNorm2(R0_2)^2);
        Q1_1 = P1_1 + bk*Q0_1; 
        Q1_2 = P1_2 + bk*Q0_2;
        R0_1 = R1_1; R0_2 = R1_2; Q0_1 = Q1_1; Q0_2 = Q1_2; X0 = X; Y0 = Y;
        k = k+1;
    end    
    Xk = X0;
end

function Y = tCTranspose(X)

% This function computes the transpose of the tensor X.

% References: P. Soto, Convergence analysis of iterative methods for computing the T-pseudoinverse
%             of complete full-rank third-order tensors based on the T-product, Results in Applied
%             Mathematics, vol. 18, p. 100372, 2023.

% Inputs: tensor X of size m x n x p

% Outputs: tensor Y of size n x m x p

    [m,n,s] = size(X);
    Y = zeros(n,m,s);
    Y(:,:,1) = (X(:,:,1))';    
    for k = 2:s
        Y(:,:,k) = (X(:,:,s-k+2))';
    end
    
end

function X = tprod3(A,B,C)

% This function computes the t-product of three third-order tensors

% Inputs: tensor A of size m x n x p
%         tensor B of size n x r x p
%         tensor C of size r x s x p

% Outputs: tensor X of size m x s x p

    Y = tprod(A,B);
    X = tprod(Y,C);
end

function C = tprod(A,B)

% This function computes the t-product of two third-order tensors

% References: P. Soto, Convergence analysis of iterative methods for computing the T-pseudoinverse
%             of complete full-rank third-order tensors based on the T-product, Results in Applied
%             Mathematics, vol. 18, p. 100372, 2023.

% Inputs: tensor A of size m x n x p
%         tensor B of size n x r x p

% Outputs: tensor C of size m x r x p

    [m1,~,p1]=size(A);
    n2=size(B,2);

    Ct=zeros(m1,n2,p1);

    At=fft(A,[],3);
    Bt=fft(B,[],3);
    halfp1=ceil((p1+1)/2);
    Ct(:,:,1:halfp1)=pagemtimes(At(:,:,1:halfp1),Bt(:,:,1:halfp1));
    idx = (halfp1+1):p1;
    mirror = p1 + 2 - idx;
    Ct(:,:,idx) = conj(Ct(:,:,mirror));
    C=real(ifft(Ct,[],3));

end

function maxi = tNorm2(A)

% This function computes the spectral norm of a third-order tensor
% under the t-product.

% References: P. Soto, Convergence analysis of iterative methods for computing the T-pseudoinverse
%             of complete full-rank third-order tensors based on the T-product, Results in Applied
%             Mathematics, vol. 18, p. 100372, 2023.

% Inputs: tensor A of size m x n x p

% Outputs: spectral norm of the tensor (maxi)

    p = size(A,3);
    maxi = -1;
    At = fft(A,[],3);
    for k = 1:p
        Sk = svd(At(:,:,k));
        SkMax = max(Sk);
        if SkMax > maxi
            maxi = SkMax;
        end
    end    

end
