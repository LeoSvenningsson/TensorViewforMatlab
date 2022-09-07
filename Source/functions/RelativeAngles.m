function [Euler] = RelativeAngles(Tensor1,Tensor2,Mode,Order1,Order2)
SymTensor = (Tensor1+transpose(Tensor1))/2; % symmetrizes the tensor
SymTensor2 = (Tensor2+transpose(Tensor2))/2; % symmetrizes the tensor
[U, D] = eig(SymTensor, 'vector'); % Rotation matrix is gained from the Eigenvectors
[U2, D2] = eig(SymTensor2, 'vector'); % Rotation matrix is gained from the Eigenvectors
Sym1 = length(D)-length(unique(round(D,7)));
Sym2 = length(D2)-length(unique(round(D2,7)));
%UnitTensor = [1 0 0; 0 1 0; 0 0 1];
% sorting the two tensors according to Order1 and Order2
if Order1=="Ascending"
    [D, ind] = sort(D,'ascend');
    U = U(:, ind);
end
if Order1=="Descending"
    [D, ind] = sort(D,'descend');
    U = U(:, ind);
end
if Order1=="Absascending"
    [~, ind] = sort(abs(D),'ascend');
    U = U(:, ind);
    D = D(ind);
end

if Order2=="Ascending"
    [D2, ind] = sort(D2,'ascend');
    U2 = U2(:, ind);
end
if Order2=="Descending"
    [D2, ind] = sort(D2,'descend');
    U2 = U2(:, ind);
end
if Order2=="Absascending"
    [~, ind] = sort(abs(D2),'ascend');
    U2 = U2(:, ind);
    D2 = D2(ind);
end



% if Sym1 == 1 && round(U(3,3),4) == 0 %Sorting only for symetrical tensors
%     U = U(:, [1 3 2]);
% end
% if Sym2 == 1 && round(U2(3,3),4) == 0 %Sorting only for symetrical tensors
%     U2 = U2(:, [1 3 2]);
% end
%%%
%Check for negative eigenvalues
PAS1(1,1)=D(1);
PAS1(2,2)=D(2);
PAS1(3,3)=D(3);
PAS2(1,1)=D2(1);
PAS2(2,2)=D2(2);
PAS2(3,3)=D2(3);

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
M=RotateTensor(Alpha,Beta,Gamma,PAS1,"AZYZ");
MF=U*PAS1*U^(-1);
if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong signa, the euler angles are incorrectly calculated
    
else
    U=-U;
end

Beta1 = acos(U2(3,3));

if U(3,3) == 1
    Alpha1 = acos(U2(1,1));
    Gamma1 = 0;
else
    Alpha1 = atan2(U2(2,3)/sin(Beta1),U2(1,3)/sin(Beta1));
    Gamma1 = atan2(U2(3,2)/sin(Beta1),-U2(3,1)/sin(Beta1));
end

Alpha = Alpha1;
Beta = Beta1;
Gamma = Gamma1;
M=RotateTensor(Alpha,Beta,Gamma,PAS2,"AZYZ");
MF=U2*PAS2*U2^(-1);

if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong sign, the euler angles are incorrectly calculated
    
else
    U2=-U2;
end
% end check negative eigen


U=U2^(-1)*U; %Sets the relative rotation matrix. Maybe an unfortunate namning consider changing

PAS(1,1)=D(1);
PAS(2,2)=D(2);
PAS(3,3)=D(3);

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
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong sign, the euler angles are incorrectly calculated
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
    if Sym1 == 1 || Sym2 == 1
        [Euler1, PASv1] = MFtoEuler(Tensor1,"AZYZ",Order1);
        [Euler2, PASv2] = MFtoEuler(Tensor2,"AZYZ",Order2);
        R1 = CreateRotationMatrix(Euler1(1),Euler1(2),Euler1(3),"AZYZ");
        R2 = CreateRotationMatrix(Euler2(1),Euler2(2),Euler2(3),"AZYZ");
        Rrel1 = R2^(-1)*R1;
        Rrel2 = R1^(-1)*R2; % not used
        
        PASv1 = round(PASv1,14);% removes floating point jittering for sensitive sqrt calculations
        PASv2 = round(PASv2,14);
        
        Arel1=round(Rrel1*PAS1*Rrel1^-1,14);
        Arel2=round(Rrel2*PAS2*Rrel2^-1,14); % not used
        
        if Sym1 == 1 && Sym2 == 1 % section for 2 symmetric tensors with a controll scheme that cheks if the tensor angle simplification is correct.
            if sum(round(Arel1(2,3),9)+round(Arel1(1,3),9)+round(Arel1(1,2),9)) == 0 %% handles the case when the relative angle is exactly [0,0,0]
                Alpha = 0; Beta = 0; Gamma = 0;
            else
                
                Alpha = 0; Beta = asin(sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZYZ");
                Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
                GammaCheck = asin(Arel1(3,2)/Mcheck(3,1)); % A trick since Mcheck(3,2) == 0 in this case
                GammarotCheck = CreateRotationMatrix(0,0,GammaCheck,"AZYZ"); % the rotaion here finds the solution on tensor A symmetry axis on a tensor B in the frame of A system.
                Mcheck=round(GammarotCheck*Mcheck*GammarotCheck^-1,3);
                
                
                if isequal(round(Arel1,3),round(Mcheck,3))
                else
                    Alpha = 0; Beta = asin(-sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                    RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZYZ");
                    Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
                    GammaCheck = asin(Arel1(3,2)/Mcheck(3,1));
                    GammarotCheck = CreateRotationMatrix(0,0,GammaCheck,"AZYZ");
                    Mcheck=round(GammarotCheck*Mcheck*GammarotCheck^-1,3);
                    if isequal(round(Arel1,3),round(Mcheck,3))
                    else
                        disp('Failed isequal check, please contact the authors for bughunting')
                    end
                end
            end
        elseif Sym1 == 1 && Sym2 == 0
            [Eulerrel1, ~] = MFtoEuler(Arel1,"AZYZ",Order1);
            Alpha = Eulerrel1(1);
            Beta = Eulerrel1(2);
            Gamma = Eulerrel1(3);
            
            RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZYZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            
            if isequal(round(Arel1,3),round(Mcheck,3))
            else
                disp('Failed isequal check, please contact the authors for bughunting')
            end
        elseif Sym1 == 0 && Sym2 == 1
            [Eulerrel1,~] = MFtoEuler(Arel2,"PZYZ",Order2);
            Alpha = Eulerrel1(1);
            Beta = Eulerrel1(2);
            Gamma = Eulerrel1(3);
            [Alpha,Beta,Gamma] = tryallanglestest(Alpha,Beta,Gamma,PASv2,PAS1,Arel1,Mode);
        end
        
    end
    if Sym1 == 2 || Sym2 == 2 % Spherical symmetry
        Alpha = 0;
        Beta = 0;
        Gamma = 0;
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
        if acos(U(3,3)) == 1
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
    
    
    if Sym1 == 1 || Sym2 == 1
        [Euler1, PASv1] = MFtoEuler(Tensor1,"AZYZ",Order1);
        [Euler2, PASv2] = MFtoEuler(Tensor2,"AZYZ",Order2);
        R1 = CreateRotationMatrix(Euler1(1),Euler1(2),Euler1(3),"AZYZ");
        R2 = CreateRotationMatrix(Euler2(1),Euler2(2),Euler2(3),"AZYZ");
        Rrel1 = R2^(-1)*R1;
        Rrel2 = R1^(-1)*R2;
        
        PASv1 = round(PASv1,12);% removes floating point jittering for sensitive sqrt calculations
        PASv2 = round(PASv2,12);
        Arel1=round(Rrel1*PAS1*Rrel1^-1,14);
        Arel2=round(Rrel2*PAS2*Rrel2^-1,14);
        
        if Sym1 == 1 && Sym2 == 1
            if sum(round(Arel1(2,3),9)+round(Arel1(1,3),9)+round(Arel1(1,2),9)) == 0 %% handles the case when the relative angle is exactly [0,0,0]
                Alpha = 0; Beta = 0; Gamma = 0;
            else
                Alpha = 0; Beta = asin(sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZYZ");
                Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
                GammaCheck = asin(Arel1(3,2)/Mcheck(3,1));
                GammarotCheck = CreateRotationMatrix(0,0,GammaCheck,"AZYZ");
                Mcheck=round(GammarotCheck*Mcheck*GammarotCheck^-1,3);
                
                %isequal(round(Arel1,3),round(Mcheck,3));
                if isequal(round(Arel1,3),round(Mcheck,3))
                else
                    Alpha = 0; Beta = asin(-sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                    RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZYZ");
                    Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
                    GammaCheck = asin(Arel1(3,2)/Mcheck(3,1));
                    GammarotCheck = CreateRotationMatrix(0,0,GammaCheck,"AZYZ");
                    Mcheck=round(GammarotCheck*Mcheck*GammarotCheck^-1,3);
                    if isequal(round(Arel1,3),round(Mcheck,3))
                    else
                        disp('Failed isequal check, please contact the authors for bughunting')
                    end
                end
            end
        elseif Sym1 == 1 && Sym2 == 0
            [Eulerrel1, ~] = MFtoEuler(Arel1,"AZYZ",Order1);
            Alpha = Eulerrel1(1);
            Beta = Eulerrel1(2);
            Gamma = Eulerrel1(3);
            
            RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZYZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            
            if isequal(round(Arel1,3),round(Mcheck,3))
            else
                disp('Failed isequal check, please contact the authors for bughunting')
            end
        elseif Sym1 == 0 && Sym2 == 1
            [Eulerrel1,~] = MFtoEuler(Arel2,"PZYZ",Order2);
            Alpha = Eulerrel1(1);
            Beta = Eulerrel1(2);
            Gamma = Eulerrel1(3);
            [Alpha,Beta,Gamma] = tryallanglestest(Alpha,Beta,Gamma,PASv2,PAS1,Arel1,Mode);
        end
        
    end
    
    if Sym1 == 2 || Sym2 == 2% Spherical symmetry
        Alpha = 0;
        Beta = 0;
        Gamma = 0;
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
        if acos(U(3,3)) == 1
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
    
    if Sym1 == 1 || Sym2 == 1
        
        [Euler1, PASv1] = MFtoEuler(Tensor1,"AZXZ",Order1);
        [Euler2, PASv2] = MFtoEuler(Tensor2,"AZXZ",Order2);
        
        R1 = CreateRotationMatrix(Euler1(1),Euler1(2),Euler1(3),"AZXZ");
        R2 = CreateRotationMatrix(Euler2(1),Euler2(2),Euler2(3),"AZXZ");
        Rrel1 = R2^(-1)*R1;
        Rrel2 = R1^(-1)*R2;
        PASv1 = round(PASv1,12);% removes floating point jittering for sensitive sqrt calculations
        PASv2 = round(PASv2,12);
        Arel1=round(Rrel1*PAS1*Rrel1^-1,14);
        Arel2=round(Rrel2*PAS2*Rrel2^-1,14);
        
        if Sym1 == 1 && Sym2 == 1
            
            if sum(round(Arel1(2,3),9)+round(Arel1(1,3),9)+round(Arel1(1,2),9)) == 0 %% handles the case when the relative angle is exactly [0,0,0]
                Alpha = 0; Beta = 0; Gamma = 0;
            else
                if  PASv1(1) == PASv1(2)
                    Alpha = pi/2; Beta = asin(sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                elseif PASv1(3) == PASv1(2)
                    Alpha = 0; Beta = pi/2; Gamma = asin(sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3))));
                end
                RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZXZ");
                Mcheck = RrelCheck*PAS1*RrelCheck^-1;
                angCheck = asin(Arel1(3,2)/Mcheck(3,1));
                angrotCheck = CreateRotationMatrix(0,0,angCheck,"AZXZ");
                Mcheck=angrotCheck*Mcheck*angrotCheck^-1;
                if isequal(round(Arel1,3),round(Mcheck,3))
                else
                    if PASv1(1) == PASv1(2)
                        Alpha = pi/2; Beta = asin(-sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                    elseif PASv1(3) == PASv1(2)
                        Alpha = 0; Beta = pi/2; Gamma = asin(-sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3))));
                    end
                    RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZXZ");
                    Mcheck = RrelCheck*PAS1*RrelCheck^-1;
                    angCheck = asin(Arel1(3,2)/Mcheck(3,1));
                    angrotCheck = CreateRotationMatrix(0,0,angCheck,"AZXZ");
                    Mcheck=angrotCheck*Mcheck*angrotCheck^-1;
                    if isequal(round(Arel1,3),round(Mcheck,3))
                    else
                        disp('Failed isequal check, please contact the authors for bughunting')
                    end
                end
            end
        elseif Sym1 == 1 && Sym2 == 0
            [Eulerrel1, ~] = MFtoEuler(Arel1,"AZXZ",Order1);
            Alpha = Eulerrel1(1);
            Beta = Eulerrel1(2);
            Gamma = Eulerrel1(3);
            
            RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZXZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            if isequal(round(Arel1,3),round(Mcheck,3))
            else
                disp('Failed isequal check at AZXZ, please contact the authors for bughunting')
            end
        elseif Sym1 == 0 && Sym2 == 1
            [Eulerrel1, ~] = MFtoEuler(Arel2,"PZXZ",Order2);
            Alpha = Eulerrel1(1);
            Beta = Eulerrel1(2);
            Gamma = Eulerrel1(3);
            [Alpha,Beta,Gamma] = tryallanglestest(Alpha,Beta,Gamma,PASv2,PAS1,Arel1,Mode);
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
    
    if Sym1 == 2 || Sym2 == 2 % setting angles for spherically symmetric tensors
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
        if acos(U(3,3)) == 1
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
    
    if Sym1 == 1 || Sym2 == 1
        [Euler1, PASv1] = MFtoEuler(Tensor1,"AZXZ",Order1);
        [Euler2, PASv2] = MFtoEuler(Tensor2,"AZXZ",Order2);
        R1 = CreateRotationMatrix(Euler1(1),Euler1(2),Euler1(3),"AZXZ");
        R2 = CreateRotationMatrix(Euler2(1),Euler2(2),Euler2(3),"AZXZ");
        Rrel1 = R2^(-1)*R1;
        Rrel2 = R1^(-1)*R2;
        
        PASv1 = round(PASv1,12);% removes floating point jittering for sensitive sqrt calculations
        PASv2 = round(PASv2,12);
        Arel1=round(Rrel1*PAS1*Rrel1^-1,14);
        Arel2=round(Rrel2*PAS2*Rrel2^-1,14);
        
        if Sym1 == 1 && Sym2 == 1
            if sum(round(Arel1(2,3),9)+round(Arel1(1,3),9)+round(Arel1(1,2),9)) == 0 %% handles the case when the relative angle is exactly [0,0,0]
                Alpha = 0; Beta = 0; Gamma = 0;
            else
                if PASv1(1) == PASv1(2) && PASv2(1) == PASv2(2)
                    Alpha = pi/2; Beta = asin(sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                elseif PASv1(3) == PASv1(2) && PASv2(3) == PASv2(2)
                    Alpha = 0; Beta = pi/2; Gamma = asin(sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3))));
                elseif PASv1(1) == PASv1(2) && PASv2(3) == PASv2(2)
                    Alpha = pi/2; Beta = asin(sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                elseif PASv1(3) == PASv1(2) && PASv2(1) == PASv2(2)
                    Alpha = 0; Beta = pi/2; Gamma = asin(sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3))));
                end
                
                RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZXZ");
                Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
                GammaCheck = asin(Arel1(3,2)/Mcheck(3,1));
                GammarotCheck = CreateRotationMatrix(0,0,GammaCheck,"AZXZ");
                Mcheck=round(GammarotCheck*Mcheck*GammarotCheck^-1,3) ;
                
                %isequal(round(Arel1,3),round(Mcheck,3));
                if isequal(round(Arel1,3),round(Mcheck,3))
                else
                    if PASv1(1) == PASv1(2) && PASv2(1) == PASv2(2)
                        Alpha = pi/2; Beta = asin(-sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                    elseif PASv1(3) == PASv1(2) && PASv2(3) == PASv2(2)
                        Alpha = 0; Beta = pi/2; Gamma = asin(-sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3))));
                    elseif PASv1(1) == PASv1(2) && PASv2(3) == PASv2(2)
                        Alpha = pi/2; Beta = asin(-sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3)))); Gamma = 0;
                    elseif PASv1(3) == PASv1(2) && PASv2(1) == PASv2(2)
                        Alpha = 0; Beta = pi/2; Gamma = asin(-sqrt((Arel1(3,3)-PASv1(3))/(PASv1(1)-PASv1(3))));
                    end
                    RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZXZ");
                    Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
                    GammaCheck = asin(Arel1(3,2)/Mcheck(3,1));
                    GammarotCheck = CreateRotationMatrix(0,0,GammaCheck,"AZXZ");
                    Mcheck=round(GammarotCheck*Mcheck*GammarotCheck^-1,3);
                    if isequal(round(Arel1,3),round(Mcheck,3))
                    else
                        disp('Failed isequal check, please contact the authors for bughunting')
                    end
                end
            end
        elseif Sym1 == 1 && Sym2 == 0
            [Eulerrel1, ~] = MFtoEuler(Arel1,"AZXZ",Order1);
            Alpha = Eulerrel1(1);
            Beta = Eulerrel1(2);
            Gamma = Eulerrel1(3);
            
            RrelCheck = CreateRotationMatrix(Alpha,Beta,Gamma,"AZXZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            if isequal(round(Arel1,3),round(Mcheck,3))
            else
                disp('Failed isequal check, please contact the authors for bughunting')
            end
        elseif Sym1 == 0 && Sym2 == 1
            [Eulerrel1, ~] = MFtoEuler(Arel2,"PZXZ",Order2);
            Alpha = Eulerrel1(1);
            Beta = Eulerrel1(2);
            Gamma = Eulerrel1(3);
            [Alpha,Beta,Gamma] = tryallanglestest(Alpha,Beta,Gamma,PASv2,PAS1,Arel1,Mode);
        end
        
    end
    
    if Sym1 == 2 || Sym2 == 2% setting angles for spherically symmetric tensors
        Alpha = 0;
        Beta = 0;
        Gamma = 0;
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
BAfinalCheck = round(U*PAS1*U^-1,14);
MfinalCheck = round(RfinalCheck*PAS1*RfinalCheck^-1,14);
if Sym1 == 0 && Sym2 == 0
    isequal(round(BAfinalCheck,3),round(MfinalCheck,3))
    
    if isequal(round(BAfinalCheck,3),round(MfinalCheck,3))
    else
        disp('Failed isequal check at RelativeRngles, please contact the authors for bughunting')
    end
end
end