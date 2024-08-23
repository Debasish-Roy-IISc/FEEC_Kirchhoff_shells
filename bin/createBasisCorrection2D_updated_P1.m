function ori=createBasisCorrection2D_updated_P1(wEdge, elType)
switch elType
    case 'P1M'
        ori=diag(wEdge);
    case 'P1'
        ori=zeros(6);
        ori(1:3,1:3)=diag(wEdge);
        ori(4:6,4:6)=diag(wEdge);

end

end