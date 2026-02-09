%%%%%%% Numerical Experiment 1 (Tables) %%%%%%%

%%%% Reference %%%%
%   Paper   = "An iterative method for computing the tensor 
%              pseudoinverse of third-order tensors using 
%              adaptable tensor-tensor products"
%   Authors = Soto-Quiros, Pablo (jusoto@tec.ac.cr) 
%             Valverde-Sanchez, Samuel (savalverde@itcr.ac.cr)

clc; clear; close all

dims=10:10:100;

timeMatrix=zeros(length(dims),6);

k=1; p=10; iterMax=100; tol=1e-10; s=3; caseNum=1;

for m=dims
    if caseNum==1    
        A=randn(m/2,3*m,10);
    elseif caseNum==2    
        B=randn(m,m/2,5); C=randn(m/2,m,5);
        A=mprod(B,C,'f');
    elseif caseNum==3    
        B=randn(m,m/2,3); C=randn(m/2,2*m,3);
        A=mprod(B,C,'f');
    end

    % Proposed Method (t-product)
    X0_t=(1/normFrob3d(A)^2)*mtranspose(A,'t');
    tic;  Xk_t=newMethod(A,X0_t,s,iterMax,tol,'t'); timeMatrix(k,1)=toc;

    % Proposed Method (c-product)
    X0_c=(1/normFrob3d(A)^2)*mtranspose(A,'c');
    tic;  Xk_c=newMethod(A,X0_c,s,iterMax,tol,'c'); timeMatrix(k,2)=toc;

    % Proposed Method (f-product)
    X0_f=(1/normFrob3d(A)^2)*mtranspose(A,'f');
    tic;  Xk_f=newMethod(A,X0_f,s,iterMax,tol,'f'); timeMatrix(k,3)=toc;

    % Karmakar Method
    tic;  Xk_k = eitpMethod(A,X0_t,5,iterMax,tol); timeMatrix(k,4)=toc;
    
    % Behera Method
    tic;  Xk_be = mhpi19Method(A,X0_t,iterMax, tol,'t'); timeMatrix(k,5)=toc;

    % Baohua Method
    X0_ba=zeros(size(A,2),size(A,1),size(A,3));
    tic;  Xk_ba = cgtpMethod(A,X0_ba,iterMax, tol); timeMatrix(k,6)=toc;
    %timeMatrix(k,6)=0;

    k=k+1;
end


methods = {'Alg. 3 (t-product)','Alg. 3 (c-product)', 'Alg. 3 (f-product)', 'EITP Method', 'MHPI19 Method', 'CGTP Method'};
T = array2table(timeMatrix, ...
    'RowNames', string(dims), ...
    'VariableNames', methods);
disp(T)