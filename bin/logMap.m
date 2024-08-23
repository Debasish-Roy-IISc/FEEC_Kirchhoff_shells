function W=logMap(R)
    A=0.5*(R-R');
    theta=(0.5*trace(A'*A))^0.5;
    W=zeros(3);
    if theta>1e-15
        W=(asin(theta)/theta)*A;
    end
end