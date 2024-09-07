function [Z,W] = loglogfit(X,Y,I)

    if(any(X)<=0 || any(Y)<=0)
        error("log(X)和log(Y)需要X和Y>0");
    end
    if(~isvector(X) || ~isvector(Y))
        error("X和Y应当是向量");
    end
    if(size(X)~=size(Y))
        error("X和Y应当拥有相同的size");
    end
    if(size(I)~=size(X))
        error("I应当和X拥有相同的size");
    end

    % 选中用于拟合的数据点
    validX = X(I);
    validY = Y(I);


    Coeff = polyfit(log(validX),log(validY),1);
    Z = polyval(Coeff,log(X));
    W = Z-Y;

    if(nargout==0)
        figure;
        loglog(X,Y);
        hold on;
        loglog(X,exp(Z));
    end

end

