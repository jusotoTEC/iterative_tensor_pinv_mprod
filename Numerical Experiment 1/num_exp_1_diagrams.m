%%%%%%% Numerical Experiment 1 (Diagrams) %%%%%%%

%%%% Reference %%%%
%   Paper   = "An iterative method for computing the tensor 
%              pseudoinverse of third-order tensors using 
%              adaptable tensor-tensor products"
%   Authors = Soto-Quiros, Pablo (jusoto@tec.ac.cr) 
%             Valverde-Sanchez, Samuel (savalverde@itcr.ac.cr)

clc; clear; close all

dims=100:100:1000;

timeMatrix=zeros(length(dims),5);
errorMatrix=zeros(length(dims),5);

k=1; iterMax=100; tol=1e-10; s=3; caseNum=1;
for m=dims
    disp(m)
    % Random Tensor

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
    tic;  [Xk_t,er_t]=newMethod(A,X0_t,s,iterMax,tol,'t'); timeMatrix(k,1)=toc;
    errorMatrix(k,1)=er_t;

    % Proposed Method (c-product)
    X0_c=(1/normFrob3d(A)^2)*mtranspose(A,'c');
    tic;  [Xk_c,er_c]=newMethod(A,X0_c,s,iterMax,tol,'c'); timeMatrix(k,2)=toc;
    errorMatrix(k,2)=er_c;

    % Proposed Method (facewise product)
    X0_f=(1/normFrob3d(A)^2)*mtranspose(A,'f');
    tic;  [Xk_f,er_f]=newMethod(A,X0_f,s,iterMax,tol,'f'); timeMatrix(k,3)=toc;
    errorMatrix(k,3)=er_f;

    % Karmakar Method
    tic;  [Xk_k,er_k] = eitpMethod(A,X0_t,5,iterMax,tol); timeMatrix(k,4)=toc;
    errorMatrix(k,4)=er_k;
    
    % Behera Method
    tic;  [Xk_b,er_b] = mhpi19Method(A,X0_t,iterMax, tol,'t'); timeMatrix(k,5)=toc;
    errorMatrix(k,5)=er_b;

    k=k+1;
end

% Plot dimension vrs time
figure
hold on
plot(dims,timeMatrix(:,1))
plot(dims,timeMatrix(:,2))
plot(dims,timeMatrix(:,3))
plot(dims,timeMatrix(:,4))
plot(dims,timeMatrix(:,5))
grid on
xlabel('Dimension (m)')
ylabel('Time (s)')
legend('Alg. 3 (t-product)','Alg. 3 (c-product)', 'Alg. 3 (f-product)', 'EITP Method', 'MHPI19 Method')
set(gca,'FontSize',18)
axis square
exportgraphics(gca,['img_time_c',num2str(caseNum),'.svg'],'ContentType','vector')

% Plot dimension vrs error
figure
hold on
loglog(dims,errorMatrix(:,1))
loglog(dims,errorMatrix(:,2))
loglog(dims,errorMatrix(:,3))
loglog(dims,errorMatrix(:,4))
loglog(dims,errorMatrix(:,5))
grid on
xlabel('Dimension (m)')
ylabel('Error')
legend('Alg. 3 (t-product)','Alg. 3 (c-product)', 'Alg. 3 (f-product)', 'EITP Method', 'MHPI19 Method')
set(gca,'FontSize',18)
axis square
exportgraphics(gca,['img_error_c',num2str(caseNum),'.svg'],'ContentType','vector')

%Percente Difference
percDiff_eitp=100*(timeMatrix(:,4)-timeMatrix(:,3))./timeMatrix(:,4);
percDiff_mhpi19=100*(timeMatrix(:,5)-timeMatrix(:,3))./timeMatrix(:,5);
figure
hold on
loglog(dims,percDiff_eitp)
loglog(dims,percDiff_mhpi19)
grid on
xlabel('Dimension (m)')
ylabel('Percent Difference')
legend('EITP Method', 'MHPI19 Method')
set(gca,'FontSize',18)
axis square
exportgraphics(gca,['img_percDiff_c',num2str(caseNum),'.svg'],'ContentType','vector')