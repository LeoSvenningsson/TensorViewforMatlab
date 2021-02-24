function  [X,Y,Z,AlphaMap] = CreateTensor(Tensor,atomcoord,CSAref,OvaloidEllipsoid,ShieldingShift,Transparency,TensorScale)

SymTensor = (Tensor+transpose(Tensor))/2; % symmetrizes the tensor

[U, D] = eig(SymTensor, 'vector'); % Rotation matrix is gained from the Eigenvectors
[D, ind] = sort(D,'ascend');
U = U(:, ind);


CShiftA = CSAref - D; % Converting chemicals shielding to shift % The app is always in "shielding" mode


Sym1 = length(D)-length(unique(round(D,7))); % checks for symetrical tensor properties

PAS(1,1)=D(1);
PAS(2,2)=D(2);
PAS(3,3)=D(3);



%%% saving the Principal Axis System shielding tensor elements in a readable variable
sxx = D(1);
syy = D(2);
szz = D(3);
%%%

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

%%%%

if Sym1 == 2 % setting angles for spherically symmetric tensors
    Alpha = 0;
    Beta = 0;
    Gamma = 0;
end

%%%%



th = (0:180)*pi/180; % resolution for the tensor surface
ph = (0:2:360)*pi/180; % resolution for the tensor surface
X = zeros(length(th),length(ph));
Y = zeros(length(th),length(ph));
Z = zeros(length(th),length(ph));

if OvaloidEllipsoid == "ovaloid" % if ovaloid
    AlphaMap = zeros(length(th),length(ph)) + Transparency; % make a map for transparency to lower transparency for negative CSA
    ORadi = zeros(length(th),length(ph));
    if ShieldingShift == 0 % if shielding
        OTensorScale=TensorScale/sum(abs(D)/3);
        for i = 1:length(th)
            for j = 1:length(ph)
                ORadi(i,j) = sxx*(sin(Gamma)*sin(Alpha - ph(j))*sin(th(i)) + cos(Gamma)*(cos(th(i))*sin(Beta) - cos(Beta)*cos(Alpha - ph(j))*sin(th(i)))).^2 + syy*(cos(th(i))*sin(Beta)*sin(Gamma) - (cos(Beta)*cos(Alpha - ph(j))*sin(Gamma) + cos(Gamma)*sin(Alpha - ph(j)))*sin(th(i))).^2 + szz*(cos(Beta)*cos(th(i)) + cos(Alpha - ph(j))*sin(Beta)*sin(th(i))).^2;
                X(i,j) = ORadi(i,j)*sin(th(i))*cos(ph(j))*OTensorScale + atomcoord(1);
                Y(i,j) = ORadi(i,j)*sin(th(i))*sin(ph(j))*OTensorScale + atomcoord(2);
                Z(i,j) = ORadi(i,j)*cos(th(i))*OTensorScale + atomcoord(3);
                CSAneg = find(ORadi<0);
                AlphaMap(CSAneg) = Transparency - 0.2; % 0.2 is the transparency shift for negative CSA
            end
        end
    elseif ShieldingShift == 1 % if shift, this part is only used when running the script version of TensorView for Matlab
        OTensorScale = TensorScale/sum(abs(CShiftA)/3);
        for i = 1:length(th)
            for j = 1:length(ph)
                ORadi(i,j) = CShiftA(1)*(sin(Gamma)*sin(Alpha - ph(j))*sin(th(i)) + cos(Gamma)*(cos(th(i))*sin(Beta) - cos(Beta)*cos(Alpha - ph(j))*sin(th(i)))).^2 + CShiftA(2)*(cos(th(i))*sin(Beta)*sin(Gamma) - (cos(Beta)*cos(Alpha - ph(j))*sin(Gamma) + cos(Gamma)*sin(Alpha - ph(j)))*sin(th(i))).^2 + CShiftA(3)*(cos(Beta)*cos(th(i)) + cos(Alpha - ph(j))*sin(Beta)*sin(th(i))).^2;
                X(i,j) = ORadi(i,j)*sin(th(i))*cos(ph(j))*OTensorScale + atomcoord(1);
                Y(i,j) = ORadi(i,j)*sin(th(i))*sin(ph(j))*OTensorScale + atomcoord(2);
                Z(i,j) = ORadi(i,j)*cos(th(i))*OTensorScale + atomcoord(3);
                CSAneg = find(ORadi<0);
                AlphaMap(CSAneg) = Transparency - 0.2; % 0.2 is the transparency shift for negative CSA
            end
        end
    else
        error("ShieldingShift must be either 0 or 1")
    end
    
elseif OvaloidEllipsoid == "ellipsoid" % if ellipsoid
    if ShieldingShift == 0 % if shielding
        ERadi = zeros(length(th),length(ph));
        ETensorScale = TensorScale/sum(abs(D)/3);
        for i = 1:length(th)
            for j = 1:length(ph)
                ERadi(i,j) = abs(sxx*syy*szz)/(sqrt(sxx^2*syy^2*(cos(th(i))*cos(Beta) + cos(ph(j) - Alpha)*sin(th(i))*sin(Beta))^2 + syy^2*szz^2*(cos(ph(j) - Alpha)*cos(Beta)*cos(Gamma)*sin(th(i)) - cos(th(i))*cos(Gamma)*sin(Beta) + sin(th(i))*sin(ph(j) - (Alpha))*sin(Gamma))^2 + sxx^2*szz^2*(cos(Gamma)*sin(th(i))*sin(ph(j) - Alpha) + (-cos(ph(j) - Alpha)*cos(Beta)*sin(th(i)) + cos(th(i))*sin(Beta))*sin(Gamma))^2));
                X(i,j) = ERadi(i,j)*sin(th(i))*cos(ph(j))*ETensorScale + atomcoord(1);
                Y(i,j) = ERadi(i,j)*sin(th(i))*sin(ph(j))*ETensorScale + atomcoord(2);
                Z(i,j) = ERadi(i,j)*cos(th(i))*ETensorScale + atomcoord(3);
            end
        end
        
        AlphaMap = zeros(length(X(:,1)),length(X(1,:))) + Transparency;
    elseif ShieldingShift == 1 % if shift
        ERadi = zeros(length(th),length(ph));
        ETensorScale = TensorScale/sum(abs(CShiftA)/3);
        for i = 1:length(th)
            for j = 1:length(ph)
                ERadi(i,j) = abs(CShiftA(1)*CShiftA(2)*CShiftA(3))/(sqrt(CShiftA(1)^2*CShiftA(2)^2*(cos(th(i))*cos(Beta) + cos(ph(j) - Alpha)*sin(th(i))*sin(Beta))^2 + CShiftA(2)^2*CShiftA(3)^2*(cos(ph(j) - Alpha)*cos(Beta)*cos(Gamma)*sin(th(i)) - cos(th(i))*cos(Gamma)*sin(Beta) + sin(th(i))*sin(ph(j) - (Alpha))*sin(Gamma))^2 + CShiftA(1)^2*CShiftA(3)^2*(cos(Gamma)*sin(th(i))*sin(ph(j) - Alpha) + (-cos(ph(j) - Alpha)*cos(Beta)*sin(th(i)) + cos(th(i))*sin(Beta))*sin(Gamma))^2));
                X(i,j) = ERadi(i,j)*sin(th(i))*cos(ph(j))*ETensorScale + atomcoord(1);
                Y(i,j) = ERadi(i,j)*sin(th(i))*sin(ph(j))*ETensorScale + atomcoord(2);
                Z(i,j) = ERadi(i,j)*cos(th(i))*ETensorScale + atomcoord(3);
            end
        end
        AlphaMap = zeros(length(X(:,1)),length(X(1,:))) + Transparency;
    else
        error("ShieldingShift must be either 0 or 1")
    end
else
    error("OvaloidEllipsoid input is not correct, use ovaloid or ellipsoid")
end


end