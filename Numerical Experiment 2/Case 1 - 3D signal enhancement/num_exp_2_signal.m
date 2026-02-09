%%%%%%% Numerical Experiment 2 (3D signal enhancement) %%%%%%%

%%%% Reference %%%%
%   Paper   = "An iterative method for computing the tensor 
%              pseudoinverse of third-order tensors using 
%              adaptable tensor-tensor products"
%   Authors = Soto-Quiros, Pablo (jusoto@tec.ac.cr) 
%             Valverde-Sanchez, Samuel (savalverde@itcr.ac.cr)

    clc; clear; close all

    a=-2; b=2;
    h1=0.1;
    h2=0.025;
    z = 0:h2:pi;

    [A,B]=tensor_Gabor(a,b,h1,h2);

    %%% Signal Reconstruction with Algorithm 2 (exact pinv) %%%
    % Compute tensor pseudoinverse of B 
    opt='f';
    Bpinv=pinvTensor(B); 
    % Construction of noisy image
    F1=mprod(A,Bpinv,opt);
    A1=mprod(F1,B,opt);

    %%% Signal Reconstruction with Algorithm 3 (f-product) %%%
    iterMax=100; tol=1e-8; s=3;
    X0=(1/normFrob3d(B)^2)*mtranspose(B,opt);
    BpinvAprox=newMethod(B,X0,s,iterMax,tol,opt); 
    % Construction of noisy image
    F2=mprod(A,BpinvAprox,opt);
    A2=mprod(F1,B,opt);

    er=normFrob3d(A1-A2);
    display(['Error given by = ', num2str(er)])

    % Generate Video
    video = VideoWriter('enhanced_signal.mp4', 'MPEG-4');
    video.FrameRate = 10;
    open(video);
    fig = figure;

    for k=1:size(A,3)
        subplot(1,4,1)
        surf(a:h1:b,a:h1:b,A(:,:,k))
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex')
        zlabel(['$f(x,y,', num2str(z(k)), ')$'], 'Interpreter', 'latex');
        title('Source Signal $f$', 'Interpreter', 'latex')
        subplot(1,4,2)
        surf(a:h1:b,a:h1:b,B(:,:,k))
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex')
        zlabel(['$h(x,y,', num2str(z(k)), ')$'], 'Interpreter', 'latex');
        title('Noisy Signal $h$', 'Interpreter', 'latex')
        subplot(1,4,3)
        surf(a:h1:b,a:h1:b,A1(:,:,k))
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex')
        zlabel(['$\widehat{f}(x,y,', num2str(z(k)), ')$'], 'Interpreter', 'latex');
        title('Reconstructed Signal $\widehat{f}$ (Alg. 2)', 'Interpreter', 'latex')
        subplot(1,4,4)
        surf(a:h1:b,a:h1:b,A2(:,:,k))
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex')
        zlabel(['$\widetilde{f}(x,y,', num2str(z(k)), ')$'], 'Interpreter', 'latex');
        title('Reconstructed Signal $\widetilde{f}$ (Alg. 3)', 'Interpreter', 'latex')
        set(gcf,'Units','normalized','Position',[0 0.25 1 0.5])
        colormap turbo

        drawnow;  % Actualiza la figura
        frame = getframe(fig);  % Captura el frame
        writeVideo(video, frame);  % Escribe el frame en el video
    end
    close(video);


