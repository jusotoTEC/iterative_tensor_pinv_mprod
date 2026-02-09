function [Xk,er] = mhpi19Method(A,X0,iterMax, tol,opt)

    [m,~,p]=size(A);
    I=meye(m,p,opt);
    Xk=X0;
    for k=1:iterMax
        Rk=I-mprod(A,Xk,opt);
        Rk2=mprod(Rk,Rk,opt);
        Uj=(7/8)*Rk+mprod(Rk2,0.5*Rk+Rk2,opt);
        Vj=(11/16)*I-(9/8)*Rk+(3/4)*Rk2+Uj;
        Aux=I+(51/128)*Rk+(39/32)*Rk2+mprod(Uj,Vj,opt);
        Xk=mprod(Xk,Aux);
        er=normFrob3d(mprod(mprod(A,Xk,opt),A,opt)-A);
        if er < tol
            break;
        end
        
    end

end