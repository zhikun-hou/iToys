function [Z] = loglogfit(X,Y)

    if(any(X)<=0 || any(Y)<=0)
        error("log(X)和log(Y)需要X和Y>0");
    end
    if(~isvector(X) || ~isvector(Y))
        error("X和Y应当是向量");
    end


    Coeff = polyfit(log(X),log(Y),1);
    Z = polyval(Coeff,log(X));

    if(nargout==0)
        figure;
        loglog(X,Y);
        hold on;
        loglog(X,exp(Z));
    end

end

