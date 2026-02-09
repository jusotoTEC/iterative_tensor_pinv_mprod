%%%%%%% Numerical Experiment 3 (Video deblurring) %%%%%%%

%%%% Reference %%%%
%   Paper   = "An iterative method for computing the tensor 
%              pseudoinverse of third-order tensors using 
%              adaptable tensor-tensor products"
%   Authors = Soto-Quiros, Pablo (jusoto@tec.ac.cr) 
%             Valverde-Sanchez, Samuel (savalverde@itcr.ac.cr)

    clc; clear; close all
    
    
    [A,fps] = video2tensor('source_video.avi');
    [B,~] = video2tensor('noisy_video.avi');     
    
    %%% Video Reconstruction with Algorithm 2 (exact pinv) %%%
    % Compute tensor pseudoinverse of B 
    opt='f';
    Bpinv=pinvTensor(B); 
    % Construction of noisy image
    F1=mprod(A,Bpinv,opt);
    A1=mprod(F1,B,opt);

    %%% Video Reconstruction with Algorithm 3 (f-product) %%%
    iterMax=100; tol=1e-8; s=3;
    X0=(1/normFrob3d(B)^2)*mtranspose(B,opt);
    BpinvAprox=newMethod(B,X0,s,iterMax,tol,opt); 
    % Construction of noisy image
    F2=mprod(A,BpinvAprox,opt);
    A2=mprod(F1,B,opt);

    er=normFrob3d(A1-A2);
    display(['Error given by = ', num2str(er)])

     % Video
    video = VideoWriter('num_exp_2_video.mp4', 'MPEG-4');
    video.FrameRate = fps;
    open(video);
    fig = figure;

    for k=1:size(A,3)
        subplot(1,4,1)
        imshow(A(:,:,k))
        title('Source Video')
        subplot(1,4,2)
        imshow(B(:,:,k))
        title('Noisy Video')
        subplot(1,4,3)
        imshow(A1(:,:,k))
        title('Reconstructed Video (Alg. 2)')
        subplot(1,4,4)
        imshow(A2(:,:,k))
        title('Reconstructed Video (Alg. 3)')
        set(gcf,'Units','normalized','Position',[0 0.25 1 0.5])


        drawnow;  % Actualiza la figura
        frame = getframe(fig);  % Captura el frame
        writeVideo(video, frame);  % Escribe el frame en el video
    end
    close(video);

