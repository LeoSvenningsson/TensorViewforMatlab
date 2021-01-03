function [Euler,PASv] = MFtoEuler(Tensor,Mode,Order)
Tensor=round(Tensor,5);
SymTensor = (Tensor+transpose(Tensor))/2; % symmetrizes the tensor
[U, D] = eig(SymTensor, 'vector'); % Rotation matrix is gained from the Eigenvectors
Sym1 = length(D)-length(unique(round(D,7)));
if Order==1
    [D, ind] = sort(D,'ascend');
    U = U(:, ind);
end
if Order==2
    [D, ind] = sort(D,'descend');
    U = U(:, ind);
end

%%%

PAS(1,1)=D(1);
PAS(2,2)=D(2);
PAS(3,3)=D(3);
PASv=D;

if Mode == 1
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
    M=RotateTensor(Alpha,Beta,Gamma,PAS,0);
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
    if Gamma>pi
        Gamma = Gamma-pi;
    end
    
    %%%%
    
    if Sym1 == 2
        Alpha = 0;
        Beta = 0;
        Gamma = 0;
    end
    
    AlphaZYZactive = Alpha;
    BetaZYZactive = Beta;
    GammaZYZactive = Gamma;
    Euler = [AlphaZYZactive,BetaZYZactive,GammaZYZactive];
end
if Mode == 2
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
    M=RotateTensor(Alpha,Beta,Gamma,PAS,0);
    MF=U*PAS*U^(-1);
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong signa, the euler angles are incorrectly calculated
        
    else
        U=-U;
        Beta1 = acos(U(3,3));
        if acos(U(3,3)) == 0
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
    if AlphaZYZpassive>pi
        AlphaZYZpassive = AlphaZYZpassive-pi;
    end
    
    if Sym1 == 2 % setting angles for spherically symmetric tensors
        AlphaZYZpassive = 0;
        BetaZYZpassive = 0;
        GammaZYZpassive = 0;
    end
    
    Euler = [AlphaZYZpassive,BetaZYZpassive,GammaZYZpassive];
end
%%% ZXZ Rotation system


%%%%
if Mode == 3
    
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
    M=RotateTensor(Alpha,Beta,Gamma,PAS,2);
    MF=U*PAS*U^(-1);
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong signa, the euler angles are incorrectly calculated
        
    else
        U=-U;
        Beta1 = acos(U(3,3));
        if acos(U(3,3)) == 0
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
    if Gamma>pi
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
if Mode == 4
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
    M=RotateTensor(Alpha,Beta,Gamma,PAS,2);
    MF=U*PAS*U^(-1);
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong signa, the euler angles are incorrectly calculated
        
    else
        U=-U;
        Beta1 = acos(U(3,3));
        if acos(U(3,3)) == 0
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
    if AlphaZXZpassive>pi
        AlphaZXZpassive = AlphaZXZpassive-pi;
    end
    
    if Sym1 == 2 % setting angles for spherically symmetric tensors
        AlphaZXZpassive = 0;
        BetaZXZpassive = 0;
        GammaZXZpassive = 0;
    end
    
    Euler = [AlphaZXZpassive,BetaZXZpassive,GammaZXZpassive];
end
end