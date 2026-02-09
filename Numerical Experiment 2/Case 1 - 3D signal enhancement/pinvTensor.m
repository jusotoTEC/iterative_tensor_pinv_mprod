function B=pinvTensor(A)

    % Tensor pseudoinverse using the f-product.

    %Reference: Jin, H., Xu, S., Wang, Y., & Liu, X. (2023). The Moore–Penrose 
    %           inverse of tensors via the M-product. Computational and Applied 
    %           Mathematics, 42(6), 294.

    [m,n,p]=size(A);
    B=zeros(n,m,p);
    for i=1:p
        B(:,:,i)=pinv(A(:,:,i));
    end
end