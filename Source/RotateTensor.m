function MF = RotateTensor(Alpha,Beta,Gamma,PAS,System)


if System == "AZYZ" || System == "PZYZ" 
    rz1 =   [cos(Alpha), -sin(Alpha), 0;
        sin(Alpha), cos(Alpha), 0;
        0, 0, 1];
    
    
    ry2 =    [cos(Beta), 0, sin(Beta);
        0, 1, 0;
        -sin(Beta), 0, cos(Beta)];
    
    
    rz3 =   [cos(Gamma), -sin(Gamma), 0;
        sin(Gamma), cos(Gamma), 0;
        0, 0, 1];
    
    
    R = rz1*ry2*rz3;
    if System == "PZYZ"
        R = rz3^(-1)*ry2^(-1)*rz1^(-1);
    end
    
elseif System == "AZXZ" || System == "PZXZ"
    
    rz1 =   [cos(Alpha), -sin(Alpha), 0;
        sin(Alpha), cos(Alpha), 0;
        0, 0, 1];
    
    rx2 =    [1, 0, 0;
        0, cos(Beta), -sin(Beta);
        0, sin(Beta), cos(Beta)];
    
    rz3 =   [cos(Gamma), -sin(Gamma), 0;
        sin(Gamma), cos(Gamma), 0;
        0, 0, 1];
    
    R = rz1*rx2*rz3;
    
    if System == "PZXZ"
        R = rz3^(-1)*rx2^(-1)*rz1^(-1);
    end
    
end

MF = R*PAS*R^(-1);

end