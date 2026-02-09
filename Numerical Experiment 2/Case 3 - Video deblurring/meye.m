function I=meye(n,p,opt)

    %       Description: Computes the identity tensor under the reduced c-product
    %   Syntax Function: I = ceye(n,p)
    %            Inputs: n = positive integer
    %                    p = positive integer 
    %            Output: I = identity tensor of dimension n x n x p
    %
    %        References: C-product toolbox                    
    %                    https://github.com/jusotoTEC/c-product-toolbox
    %
    %   Code written by: Pablo Soto-Quiros (jusoto@tec.ac.cr) and
    %                    Samuel Valverde-Sanchez (savalverde@itcr.ac.cr)
    if strcmpi(opt,'t')
        I=zeros(n,n,p);
        I(:,:,1)=eye(n);
    elseif strcmpi(opt,'c')
        It=repmat(eye(n), 1, 1, p);   
        I=idct(It, [], 3);     
    elseif strcmpi(opt,'f')
        I = repmat(eye(n), 1, 1, p);
    elseif isnumeric(opt) && ismatrix(opt) && all(size(opt) >= [2 2])
        M=opt;
        It=repmat(eye(n), 1, 1, p);   
        I=mode3_product(It,M^-1);
    else
        error("The third input must be 't' (t-product), 'c' (reduced c-product), 'f' (facewise product), or a matrix M.")    
    end
end