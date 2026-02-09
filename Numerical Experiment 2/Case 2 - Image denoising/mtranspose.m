function B=mtranspose(A,opt)
    if strcmpi(opt,'t')
        [n1,n2,n3] = size(A);
        B = zeros(n2,n1,n3);
        B(:,:,1) = A(:,:,1)';
        indx=2:n3;
        mirror=n3-indx+2;
        B(:,:,indx) = permute(A(:,:,mirror),[2 1 3]);
    elseif strcmpi(opt,'c')
        At=dct(A,[],3);
        Bt=permute(At, [2 1 3]);
        B=idct(Bt,[],3);
    elseif strcmpi(opt,'f')
        B = permute(A, [2 1 3]);
    elseif isnumeric(opt) && ismatrix(opt) && all(size(opt) >= [2 2])
        M=opt;
        At=mode3_product(A,M);
        Bt=permute(At, [2 1 3]);
        B=mode3_product(Bt,M^-1);
    else
        error("The third input must be 't' (t-product), 'c' (reduced c-product), 'f' (facewise product), or a matrix M.")
    end
end