function sqEulerR =  sqminEuler(D,SymTensor,XY) % obtain square min function
 %%% Do a few random point searches to avoid local minima
 fval=1000;

if XY == "ZYZ"
    while fval>10^(-1)
    randEuler=rand(1,2).*[2*pi,pi/2];
    fval=EulerRY(randEuler,D,SymTensor);
    end
end
if XY == "ZXZ"
    while fval>10^(-1)
    randEuler=rand(1,2).*[2*pi,pi/2];
    fval=EulerRX(randEuler,D,SymTensor);
    end
end
 %%%
 
%%% Iterative min function search with the random point solution as starting values
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0,0];
x0 = randEuler;
if XY == "ZYZ"
    ub = [2*pi,pi/2];
    EulerR = @(r) EulerRY(r,D,SymTensor);
    [Euler,Fval] = fmincon(EulerR,x0,A,b,Aeq,beq,lb,ub);
    sqEulerR = [Euler 0];
end
if XY == "ZXZ"
    ub = [2*pi,pi/2];
    EulerR = @(r) EulerRX(r,D,SymTensor);
    [Euler,Fval] = fmincon(EulerR,x0,A,b,Aeq,beq,lb,ub);
    sqEulerR = [0 Euler(2) Euler(1)];
end


%%
end
%Functon used to find the Euler angles from the Principal Axis System to the Molecular Frame CSA tensor
function CSAEulerRmin = EulerRY(r,D,SymTensor)
        rza =   [cos(r(1)), -sin(r(1)), 0;
            sin(r(1)), cos(r(1)), 0;
            0, 0, 1];
        
        ryb =  [cos(r(2)), 0, sin(r(2));
            0, 1, 0;
            -sin(r(2)), 0, cos(r(2))];
        
        rzg =  [cos(0), -sin(0), 0;
            sin(0), cos(0), 0;
            0, 0, 1];

        CSAEulerR =(rza*ryb*rzg)*[D(1), 0, 0; 0, D(2), 0; 0 ,0 ,D(3)]*(rza*ryb*rzg)^(-1);
        
        CSAEulerRmin = sum(sum((CSAEulerR-SymTensor).^2)); % min square eq
end
function CSAEulerRmin = EulerRX(r,D,SymTensor)
        rza =  [cos(0), -sin(0), 0;
                sin(0), cos(0), 0;
                0, 0, 1];
        
        
        ryb = [1, 0, 0;
               0, cos(r(2)), -sin(r(2));
               0, sin(r(2)), cos(r(2))];
        
        rzg =   [cos(r(1)), -sin(r(1)), 0;
                sin(r(1)), cos(r(1)), 0;
                0, 0, 1];

        CSAEulerR =(rza*ryb*rzg)*[D(1), 0, 0; 0, D(2), 0; 0 ,0 ,D(3)]*(rza*ryb*rzg)^(-1);
        
        CSAEulerRmin = sum(sum((CSAEulerR-SymTensor).^2)); % min square eq
end