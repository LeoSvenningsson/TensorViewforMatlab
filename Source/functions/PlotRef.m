function  PlotRef(Angle,x,y,z,Scale,Brightness,Mode)
System = Mode;
Alpha = Angle(1);
Beta = Angle(2);
Gamma = Angle(3);

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

Xref = R*[1;0;0]*Scale;
Yref = R*[0;1;0]*Scale;
Zref = R*[0;0;1]*Scale;


BrR1 = min([1 (1+Brightness)]);
BrG1 = min([1 (Brightness)]);
BrB1 = min([1 (Brightness)]);

BrG1 = max([0 (Brightness)]);
BrB1 = max([0 (Brightness)]);

BrR2 = min([1 (Brightness)]);
BrG2 = min([1 (1+Brightness)]);
BrB2 = min([1 (Brightness)]);

BrR2 = max([0 (Brightness)]);
BrB2 = max([0 (Brightness)]);

BrR3 = min([1 (Brightness)]);
BrG3 = min([1 (Brightness)]);
BrB3 = min([1 (1+Brightness)]);

BrR3 = max([0 (Brightness)]);
BrG3 = max([0 (Brightness)]);

hold on
quiver3(x,y,z,Xref(1),Xref(2),Xref(3),'Color',[BrR1,BrG1,BrB1],'LineWidth',3,'AutoScale','off')

axis equal
set(gcf,'Color','w')
set(gca,'visible','off')
quiver3(x,y,z,Yref(1),Yref(2),Yref(3),'Color',[BrR2,BrG2,BrB2],'LineWidth',3,'AutoScale','off')
quiver3(x,y,z,Zref(1),Zref(2),Zref(3),'Color',[BrR3,BrG3,BrB3],'LineWidth',3,'AutoScale','off')
%hold off


end