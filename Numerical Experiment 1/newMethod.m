function [Xk,er]=newMethod(A,Xk,s,iterMax,tol,opt)

% References: Valverde-Sanchez, Samuel and Soto-Quiros, Pablo. An iterative 
%             method for computing the tensor pseudoinverse of third-order 
%             tensors using adaptable tensor-tensor products.

    [m,~,p]=size(A);
    r = 1:s;
    c = (-1).^(r-1) .* arrayfun(@(rr) nchoosek(s, rr), r);
    P=mprod(A,Xk,opt);
    I=meye(m,p,opt);
    for k=0:iterMax
        Y=c(s)*I;        
        for r=s-1:-1:1
            Y=mprod(Y,P,opt)+c(r)*I;
        end
        Xk=mprod(Xk,Y,opt);
        P=mprod(A,Xk,opt);
        Aux=mprod(P,A,opt);
        er=normFrob3d(Aux-A);
        if er < tol
            break;
        end
    end
end