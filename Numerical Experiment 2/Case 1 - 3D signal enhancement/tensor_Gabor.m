function [A,C]=tensor_Gabor(a,b,h1,h2)
    x = a:h1:b;
    y = a:h1:b;
    z = 0:h2:pi;
    s=length(z);
    A=zeros(length(x),length(y),s);
    for k=1:s
        A(:,:,k) = generar_matriz(x, y, z(k));
    end
    N=rand(size(A));
    C=A.*N;      
end

function H = generar_matriz(x,y,p)
    u = 0.25; 
    [X, Y] = meshgrid(x, y);
    H = exp(-(X.^2 + Y.^2)/2) .* cos(2*pi*u*X + p);
end