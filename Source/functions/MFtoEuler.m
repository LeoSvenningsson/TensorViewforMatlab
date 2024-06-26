function [Euler,PASv] = MFtoEuler(Tensor,Mode,Order)
SymTensor = (Tensor+transpose(Tensor))/2; % symmetrizes the tensor
[U, D] = eig(SymTensor, 'vector'); % Rotation matrix is gained from the Eigenvectors
Sym1 = length(D)-length(unique(round(D,7)));

if Order=="Ascending"
    [D, ind] = sort(D,'ascend');
    U = U(:, ind);
end
if Order=="Descending"
    [D, ind] = sort(D,'descend');
    U = U(:, ind);
end
if Order=="Absascending"
    [~, ind] = sort(abs(D),'ascend');
    U = U(:, ind);
    D = D(ind);
end

%if Sym1 == 1 && round(U(3,3),4) == 0 %Sorting only for symetrical tensors
%    U = U(:, [1 3 2]);
%end
%%%

PAS(1,1)=D(1);
PAS(2,2)=D(2);
PAS(3,3)=D(3);
PASv=D;

if Mode == "AZYZ"
    %%% ZYZ Active Rotation system
    
    Beta1 = acos(U(3,3));
    
    if U(3,3) == 1
        Alpha1 = acos(U(1,1));
        Gamma1 = 0;
    else
        Alpha1 = atan2(U(2,3)/sin(Beta1),U(1,3)/sin(Beta1));
        Gamma1 = atan2(U(3,2)/sin(Beta1),-U(3,1)/sin(Beta1));
    end
    
    Alpha = Alpha1;
    Beta = Beta1;
    Gamma = Gamma1;
    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
    MF=U*PAS*U^(-1);
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong signa, the euler angles are incorrectly calculated
        
    else
        U=-U;
        Beta1 = acos(U(3,3));
        if U(3,3) == 1
            Alpha1 = acos(U(1,1));
            Gamma1 = 0;
        else
            Alpha1 = atan2(U(2,3)/sin(Beta1),U(1,3)/sin(Beta1));
            Gamma1 = atan2(U(3,2)/sin(Beta1),-U(3,1)/sin(Beta1));
        end
        Alpha = Alpha1;
        Beta = Beta1;
        Gamma = Gamma1;
    end
    
    %%%
    if Sym1 == 2 % Spherical symmetry
        Alpha = 0;
        Beta = 0;
        Gamma = 0;
    end
    if Sym1 == 1
        if round(D(3),4) == round(D(2),4)
            %sqEul  =  sqminEuler(D,SymTensor,"ZYZ"); % Alternative calc
            %Alpha = sqEul(1);
            %Beta = sqEul(2);
            %Gamma = sqEul(3);
            SymTensor = round(SymTensor,14); % removes floating point jittering
            D = round(D,14);
            if round(SymTensor(2,3),4) == 0 && round(SymTensor(1,2),4) == 0
                if abs(asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))))) == pi/2
                    Gamma = asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                    Beta = 0;
                    Alpha = 0;
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                        Beta = 0;
                        Alpha = 0;
                    end
                    
                else
                    Gamma = asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                    Beta = asin(sqrt((SymTensor(3,3)-D(3))/(D(1)-D(3) + (D(2)-D(1))*sin(Gamma)^2)));
                    Alpha = 0;
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                        Beta = asin(sqrt((SymTensor(3,3)-D(3))/(D(1)-D(3) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                        Beta = asin(-sqrt((SymTensor(3,3)-D(3))/(D(1)-D(3) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                        Beta = asin(-sqrt((SymTensor(3,3)-D(3))/(D(1)-D(3) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                end
            else
                Gamma = asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                Beta = atan2(-SymTensor(2,3)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))),SymTensor(1,2)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))));
                Alpha = 0;
                
                M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                if isequal(round(M,3),round(SymTensor,3))
                else
                    Gamma = asin(-sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                    Beta = atan2(-SymTensor(2,3)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))),SymTensor(1,2)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))));
                    Alpha = 0;
                end
            end
        else
            Gamma = 0;
        end
    end
    Alpha = mod(Alpha,2*pi);
    Beta = mod(Beta,2*pi);
    Gamma = mod(Gamma,2*pi);
    
    if Beta>pi
        Alpha = Alpha - pi;
        Alpha = mod(Alpha,2*pi);
        Beta =  2*pi - Beta;
    end
    
    if Beta>=pi/2
        Alpha = pi + Alpha;
        Alpha = mod(Alpha,2*pi);
        Beta = pi - Beta;
        Beta = mod(Beta,2*pi);
        Gamma =  pi - Gamma;
        Gamma = mod(Gamma,2*pi);
    end
    if Gamma>=pi
        Gamma = Gamma-pi;
    end
    
    %%%%
    
    AlphaZYZactive = Alpha;
    BetaZYZactive = Beta;
    GammaZYZactive = Gamma;
    Euler = [AlphaZYZactive,BetaZYZactive,GammaZYZactive];
end
if Mode == "PZYZ"
    %%% ZYZ Passive Rotation system
    Beta1 = acos(U(3,3));
    
    if U(3,3) == 1
        Alpha1 = acos(U(1,1));
        Gamma1 = 0;
    else
        Alpha1 = atan2(U(2,3)/sin(Beta1),U(1,3)/sin(Beta1));
        Gamma1 = atan2(U(3,2)/sin(Beta1),-U(3,1)/sin(Beta1));
    end
    
    Alpha = Alpha1;
    Beta = Beta1;
    Gamma = Gamma1;
    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
    MF=U*PAS*U^(-1);
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong signa, the euler angles are incorrectly calculated
        
    else
        U=-U;
        Beta1 = acos(U(3,3));
        if U(3,3) == 1
            Alpha1 = acos(U(1,1));
            Gamma1 = 0;
        else
            Alpha1 = atan2(U(2,3)/sin(Beta1),U(1,3)/sin(Beta1));
            Gamma1 = atan2(U(3,2)/sin(Beta1),-U(3,1)/sin(Beta1));
        end
        Alpha = Alpha1;
        Beta = Beta1;
        Gamma = Gamma1;
    end
    
    if Sym1 == 2 % Spherical symmetry
        Alpha = 0;
        Beta = 0;
        Gamma = 0;
    end
    
    if Sym1 == 1
        if round(D(3),4) == round(D(2),4)
            %sqEul  =  sqminEuler(D,SymTensor,"ZYZ"); % Alternative calc
            %Alpha = sqEul(1);
            %Beta = sqEul(2);
            %Gamma = sqEul(3);
            %SymTensor;
            SymTensor = round(SymTensor,14); % removes floating point jittering
            D = round(D,14);
            if round(SymTensor(2,3),4) == 0 && round(SymTensor(1,2),4) == 0
                if abs(asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))))) == pi/2
                    
                    Gamma = asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                    Beta = 0;
                    Alpha = 0;
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                        Beta = 0;
                        Alpha = 0;
                    end
                    
                else
                    Gamma = asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                    Beta = asin(sqrt((SymTensor(3,3)-D(3))/(D(1)-D(3) + (D(2)-D(1))*sin(Gamma)^2)));
                    Alpha = 0;
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                        Beta = asin(sqrt((SymTensor(3,3)-D(3))/(D(1)-D(3) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                        Beta = asin(-sqrt((SymTensor(3,3)-D(3))/(D(1)-D(3) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                        Beta = asin(-sqrt((SymTensor(3,3)-D(3))/(D(1)-D(3) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                end
            else
                Gamma = asin(sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                Beta = atan2(-SymTensor(2,3)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))),SymTensor(1,2)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))));
                Alpha = 0;
                
                M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
                if isequal(round(M,3),round(SymTensor,3))
                else
                    Gamma = asin(-sqrt((SymTensor(2,2)-D(2))/(D(1)-D(2))));
                    Beta = atan2(-SymTensor(2,3)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))),SymTensor(1,2)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))));
                    Alpha = 0;
                end
            end
        else
            Gamma = 0;
        end
    end

    AlphaZYZpassive = mod(-Gamma,2*pi);
    BetaZYZpassive = mod(-Beta,2*pi);
    GammaZYZpassive = mod(-Alpha,2*pi);
    
    if BetaZYZpassive>pi
        BetaZYZpassive =  2*pi - BetaZYZpassive;
        GammaZYZpassive = GammaZYZpassive - pi;
        GammaZYZpassive = mod(GammaZYZpassive,2*pi);
    end
    if BetaZYZpassive>=pi/2
        AlphaZYZpassive = - (AlphaZYZpassive - pi);
        AlphaZYZpassive = mod(AlphaZYZpassive,2*pi);
        BetaZYZpassive = - (BetaZYZpassive - pi);
        BetaZYZpassive = mod(BetaZYZpassive,2*pi);
        GammaZYZpassive = GammaZYZpassive + pi;
        GammaZYZpassive = mod(GammaZYZpassive,2*pi);
    end
    if AlphaZYZpassive>=pi
        AlphaZYZpassive = AlphaZYZpassive-pi;
    end
    
    Euler = [AlphaZYZpassive,BetaZYZpassive,GammaZYZpassive];
end
%%% ZXZ Rotation system


%%%%
if Mode == "AZXZ" % Active ZXZ
    
    Beta1 = acos(U(3,3));
    
    if U(3,3) == 1
        Alpha1 = acos(U(1,1));
        Gamma1 = 0;
    else
        Alpha1 = atan2(U(1,3)/sin(Beta1),-U(2,3)/sin(Beta1));
        Gamma1 = atan2(U(3,1)/sin(Beta1),U(3,2)/sin(Beta1));
    end
    
    Alpha = Alpha1;
    Beta = Beta1;
    Gamma = Gamma1;
    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
    MF=U*PAS*U^(-1);
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong signa, the euler angles are incorrectly calculated
        
    else
        U=-U;
        Beta1 = acos(U(3,3));
        if U(3,3) == 1
            Alpha1 = acos(U(1,1));
            Gamma1 = 0;
        else
            Alpha1 = atan2(U(1,3)/sin(Beta1),-U(2,3)/sin(Beta1));
            Gamma1 = atan2(U(3,1)/sin(Beta1),U(3,2)/sin(Beta1));
        end
        Alpha = Alpha1;
        Beta = Beta1;
        Gamma = Gamma1;
    end
    
    if Sym1 == 1
        if round(D(3),4) == round(D(2),4)
            %sqEul  =  sqminEuler(D,SymTensor,"ZXZ"); % Alternative calc
            %Alpha = sqEul(1)
            %Beta = sqEul(2)
            %Gamma = sqEul(3)
            SymTensor = round(SymTensor,14); % removes floating point jittering for sensitive sqrt operations
            D = round(D,14); 
            if round(SymTensor(3,1),4) == 0 && round(SymTensor(2,1),4) == 0
                if asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1)))) == 0 || asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1)))) == pi
                    Gamma = asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                    Beta = 0;
                    Alpha = 0;
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                        Beta = 0;
                        Alpha = 0;
                    end
                    
                else
                    Gamma = asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                    Beta = asin(sqrt((SymTensor(2,2)-D(2) - (D(1)-D(2))*sin(Gamma)^2)/(D(3)-D(2) + (D(2)-D(1))*sin(Gamma)^2)));
                    Alpha = 0;
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                    if isequal(round(M,3),round(SymTensor,3)) % Perhaps i dont need to check all permutations of sqrt signs, but its in the just in case.
                    else
                        Gamma = asin(-sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                        Beta = asin(sqrt((SymTensor(2,2)-D(2) - (D(1)-D(2))*sin(Gamma)^2)/(D(3)-D(2) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                        Beta = asin(-sqrt((SymTensor(2,2)-D(2) - (D(1)-D(2))*sin(Gamma)^2)/(D(3)-D(2) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                        Beta = asin(-sqrt((SymTensor(2,2)-D(2) - (D(1)-D(2))*sin(Gamma)^2)/(D(3)-D(2) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                end
            else
                Gamma = asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                Beta = atan2(SymTensor(3,1)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))),SymTensor(2,1)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))));
                Alpha = 0;
                
                M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                
                if isequal(round(M,3),round(SymTensor,3))
                else
                    Gamma = asin(-sqrt((SymTensor(1,1)-D(2))/(D(2)-D(1))));
                    Beta = atan2(SymTensor(3,1)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))),SymTensor(2,1)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))));
                    Alpha = 0;
                end
            end
        else
            Gamma = 0;
        end
    end
    
    Alpha = mod(Alpha,2*pi);
    Beta = mod(Beta,2*pi);
    Gamma = mod(Gamma,2*pi);
    
    if Beta>pi
        Alpha = Alpha - pi;
        Alpha = mod(Alpha,2*pi);
        Beta =  2*pi - Beta;
    end
    
    if Beta>=pi/2
        Alpha = pi + Alpha;
        Alpha = mod(Alpha,2*pi);
        Beta = pi - Beta;
        Beta = mod(Beta,2*pi);
        Gamma =  pi - Gamma;
        Gamma = mod(Gamma,2*pi);
    end
    if Gamma>=pi
        Gamma = Gamma-pi;
    end
    
    if Sym1 == 2 % setting angles for spherically symmetric tensors
        Alpha = 0;
        Beta = 0;
        Gamma = 0;
    end
    
    
    
    AlphaZXZactive = Alpha;
    BetaZXZactive = Beta;
    GammaZXZactive = Gamma;
    
    Euler = [AlphaZXZactive,BetaZXZactive,GammaZXZactive];
end
if Mode == "PZXZ" % Passive ZYZ
    Beta1 = acos(U(3,3));
    
    if U(3,3) == 1
        Alpha1 = acos(U(1,1));
        Gamma1 = 0;
    else
        Alpha1 = atan2(U(1,3)/sin(Beta1),-U(2,3)/sin(Beta1));
        Gamma1 = atan2(U(3,1)/sin(Beta1),U(3,2)/sin(Beta1));
    end
    
    Alpha = Alpha1;
    Beta = Beta1;
    Gamma = Gamma1;
    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
    MF=U*PAS*U^(-1);
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong signa, the euler angles are incorrectly calculated
        
    else
        U=-U;
        Beta1 = acos(U(3,3));
        if U(3,3) == 1
            Alpha1 = acos(U(1,1));
            Gamma1 = 0;
        else
            Alpha1 = atan2(U(1,3)/sin(Beta1),-U(2,3)/sin(Beta1));
            Gamma1 = atan2(U(3,1)/sin(Beta1),U(3,2)/sin(Beta1));
        end
        Alpha = Alpha1;
        Beta = Beta1;
        Gamma = Gamma1;
    end
    
    if Sym1 == 2 % setting angles for spherically symmetric tensors
        Alpha = 0;
        Beta = 0;
        Gamma = 0;
    end
    if Sym1 == 1
        if round(D(3),4) == round(D(2),4)
            %sqEul  =  sqminEuler(D,SymTensor,"ZXZ");
            %Alpha = sqEul(1);
            %Beta = sqEul(2);
            %Gamma = sqEul(3);
            SymTensor = round(SymTensor,14); % removes floating point jittering for sensitive sqrt operations
            D = round(D,14); 
            if round(SymTensor(3,1),4) == 0 && round(SymTensor(2,1),4) == 0
                if asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1)))) == 0 || asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1)))) == pi
                    Gamma = asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                    Beta = 0;
                    Alpha = 0;
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                        Beta = 0;
                        Alpha = 0;
                    end
                    
                else
                    Gamma = asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                    Beta = asin(sqrt((SymTensor(2,2)-D(2) - (D(1)-D(2))*sin(Gamma)^2)/(D(3)-D(2) + (D(2)-D(1))*sin(Gamma)^2)));
                    Alpha = 0;
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                        Beta = asin(sqrt((SymTensor(2,2)-D(2) - (D(1)-D(2))*sin(Gamma)^2)/(D(3)-D(2) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                        Beta = asin(-sqrt((SymTensor(2,2)-D(2) - (D(1)-D(2))*sin(Gamma)^2)/(D(3)-D(2) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                    if isequal(round(M,3),round(SymTensor,3))
                    else
                        Gamma = asin(-sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                        Beta = asin(-sqrt((SymTensor(2,2)-D(2) - (D(1)-D(2))*sin(Gamma)^2)/(D(3)-D(2) + (D(2)-D(1))*sin(Gamma)^2)));
                        Alpha = 0;
                    end
                end
            else
                Gamma = asin(sqrt((SymTensor(1,1)-D(1))/(D(2)-D(1))));
                Beta = atan2(SymTensor(3,1)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))),SymTensor(2,1)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))));
                Alpha = 0;
                
                M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZXZ");
                if isequal(round(M,3),round(SymTensor,3))
                else
                    Gamma = asin(-sqrt((SymTensor(1,1)-D(2))/(D(2)-D(1))));
                    Beta = atan2(SymTensor(3,1)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))),SymTensor(2,1)/(sin(Gamma)*cos(Gamma)*(D(1)-D(2))));
                    Alpha = 0;
                end
            end
        else
            Gamma = 0;
        end
    end
    
    AlphaZXZpassive = mod(-Gamma,2*pi);
    BetaZXZpassive = mod(-Beta,2*pi);
    GammaZXZpassive = mod(-Alpha,2*pi);
    
    if BetaZXZpassive>pi
        BetaZXZpassive =  2*pi - BetaZXZpassive;
        GammaZXZpassive = GammaZXZpassive - pi;
        GammaZXZpassive = mod(GammaZXZpassive,2*pi);
    end
    if BetaZXZpassive>=pi/2
        AlphaZXZpassive = pi - AlphaZXZpassive;
        AlphaZXZpassive = mod(AlphaZXZpassive,2*pi);
        BetaZXZpassive = pi - BetaZXZpassive;
        BetaZXZpassive = mod(BetaZXZpassive,2*pi);
        GammaZXZpassive = GammaZXZpassive + pi;
        GammaZXZpassive = mod(GammaZXZpassive,2*pi);
    end
    if AlphaZXZpassive>=pi
        AlphaZXZpassive = AlphaZXZpassive-pi;
    end

    Euler = [AlphaZXZpassive,BetaZXZpassive,GammaZXZpassive];
end

% a Final check if the tensor angles are correct
RfinalCheck = CreateRotationMatrix(Euler(1),Euler(2),Euler(3),Mode);
UfinalCheck = round(U*PAS*U^-1,14);
MfinalCheck = round(RfinalCheck*PAS*RfinalCheck^-1,14);

if isequal(round(UfinalCheck,3),round(MfinalCheck,3))
else
    disp('Failed isequal check at MFtoEuler, please contact the authors for bughunting')
end

end
