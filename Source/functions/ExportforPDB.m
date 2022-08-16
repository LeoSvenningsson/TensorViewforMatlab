function ExportforPDB(fileName,TensorList)
fileID = fopen(fileName,'w');
for ntens = 1:length(TensorList(1,:))
    
    Tensor = zeros(3,3);
    
    Tensor(1,1) = TensorList((1),ntens);
    Tensor(1,2) = TensorList((2),ntens);
    Tensor(1,3) = TensorList((3),ntens);
    Tensor(2,1) = TensorList((4),ntens);
    Tensor(2,2) = TensorList((5),ntens);
    Tensor(2,3) = TensorList((6),ntens);
    Tensor(3,1) = TensorList((7),ntens);
    Tensor(3,2) = TensorList((8),ntens);
    Tensor(3,3) = TensorList((9),ntens);
    
    
    
    SymTensor = (Tensor+transpose(Tensor))/2; % symmetrizes the tensor
    
    [U, D] = eig(SymTensor, 'vector'); % Rotation matrix is gained from the Eigenvectors
    [D, ind] = sort(D,'ascend');
    U = U(:, ind);
    
    Sym1 = length(D)-length(unique(round(D,7))); % checks for symetrical tensor properties
    
    PAS(1,1)=D(1);
    PAS(2,2)=D(2);
    PAS(3,3)=D(3);
    
    
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
    M=RotateTensor(Alpha,Beta,Gamma,PAS,"AZYZ");
    MF=U*PAS*U^(-1);
    
    if isequal(round(M,3),round(MF,3)) % This check is to solve the cases when "eig" produces negative eigen vectors. If the eigen vectors have the wrong sign, the euler angles are incorrectly calculated
        
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
    
    
    s=1;% default from Leonard Muells python conversion
    
    abseig = abs(D);
    
    q=max(abseig);
    
    Mq=q^2;
    
    sS=s^2;
    
    DD = D.^2;
    
    diagDD = [DD(1) 0 0; 0 DD(2) 0; 0 0 DD(3)];
    
    matrixM = U*diagDD*(U^-1);
    
    matrixS=round(matrixM*sS*(10000/Mq),5);
    matrixSrow = [ntens matrixS(1,1) matrixS(2,2) matrixS(3,3) matrixS(1,2) matrixS(1,3) matrixS(2,3)]
    
    fprintf(fileID,'ANISOU    %1d  X           0   %-7.0f%-7.0f%-7.0f%-7.0f%-7.0f%-13.0fX\n',matrixSrow);
    
end
fclose(fileID);
end