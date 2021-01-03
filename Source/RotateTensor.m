function MF = RotateTensor(Alpha,Beta,Gamma,PAS,System)


if System == 0 || System == 1 % 0 = ZYZactive; 1 = ZYZpassive; 2 = ZXZactive; 3 = ZXZpassive
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
    if System == 1
        R = rz3^(-1)*ry2^(-1)*rz1^(-1);
    end
    
elseif System == 2 || System == 3
    
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
    
    if System == 3
        R = rz3^(-1)*rx2^(-1)*rz1^(-1);
    end
    
end

MF = R*PAS*R^(-1);

end