%%%%%%% Numerical Experiment 2 (Image denoising) %%%%%%%

%%%% Reference %%%%
%   Paper   = "An iterative method for computing the tensor 
%              pseudoinverse of third-order tensors using 
%              adaptable tensor-tensor products"
%   Authors = Soto-Quiros, Pablo (jusoto@tec.ac.cr) 
%             Valverde-Sanchez, Samuel (savalverde@itcr.ac.cr)

clc; clear; close all

%Source image
numExample=3;
if numExample==1
    A=im2double(imread("waterfall.png"));
elseif numExample==2
    A=im2double(imread("basilica.png"));
elseif numExample==3
    A=im2double(imread("beach.png"));
end

subplot(1,4,1)
imshow(A)
title('Source image')

%Noisy image
N=randn(size(A));
B=A+0.3*N;
subplot(1,4,2)
imshow(B)
title('Noisy image')

%%%% Reconstruction with Algorithm 3 (f-product) %%%%

% Approximate tensor pseudoinverse of B 
iterMax=100; tol=1e-8; s=3; opt='f';
X0=(1/normFrob3d(B)^2)*mtranspose(B,opt);
BpinvAprox=newMethod(B,X0,s,iterMax,tol,opt); 
% Construction of noisy image
F1=mprod(A,BpinvAprox,opt);
A_rec_1=mprod(F1,B,opt);
subplot(1,4,3)
imshow(A_rec_1)
title('Denoising image (Algorithm 2)')

%%%% Reconstruction with exact formula of pinv (f-product) %%%%

% Compute tensor pseudoinverse of B 
Bpinv=pinvTensor(B); 
% Construction of noisy image
F2=mprod(A,Bpinv,opt);
A_rec_2=mprod(F2,B,opt);
subplot(1,4,4)
imshow(A_rec_2)
title('Denoising image (Algorithm 3)')


% Compute color SSIM between A_rec_1 y A_rec_2
ssimR=ssim(A_rec_1(:,:,1),A_rec_2(:,:,1));
ssimG=ssim(A_rec_1(:,:,2),A_rec_2(:,:,2));
ssimB=ssim(A_rec_1(:,:,3),A_rec_2(:,:,3));

ssimRGB=(1/3)*(ssimR+ssimG+ssimB);
disp(['Color SSIM is ', num2str(ssimRGB)])



