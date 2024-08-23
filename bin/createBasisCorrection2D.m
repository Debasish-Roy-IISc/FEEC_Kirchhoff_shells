function ori=createBasisCorrection2D(wEdge, elType)
switch elType
    case 'P1M'
        ori=diag(wEdge);
    case 'P1'
        ori=eye(6);
        for i=1:3
            if wEdge(i)==-1
                ori(i,i)=0;
                ori(3+i,3+i)=0;
                ori(i,3+i)=1;
                ori(3+i,i)=1;
            end
        end
end

end