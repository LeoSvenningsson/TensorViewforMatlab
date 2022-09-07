function [AlphaOut,BetaOut,GammaOut] = tryallanglestest(Alpha,Beta,Gamma,PASv2,PAS1,Arel1,Mode)
% the function tests all equivalent angles of the symmetric tensor
% "tensor 2". We go through the 4 passive angles here since we use the trick
% to calculate tensor 1 relative angles in the frame of tensor 2 form the
% passive angles of tensor 2 in the frame of tensor 1. It might look
% strange to do this corner case relative rotation in this way, but since
% we use the pipline that handles symmetric tensors we obtain the solution with one relative
% angle at zero.
if Mode == "AZYZ" || Mode == "PZYZ"
    if PASv2(1) == PASv2(2)
        AlphaOut = Alpha; BetaOut = Beta; GammaOut = Gamma;
        RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZYZ");
        Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
        component1 = Arel1(3,1) + Arel1(3,2)*Mcheck(3,2)/Mcheck(3,1);
        component2 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
        component3 = Arel1(3,2) - Arel1(3,1)*Mcheck(3,2)/Mcheck(3,1);
        component4 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
        %SymrotangCheckcos = acos(component1/component2) %a Z rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
        %SymrotangChecksin = asin(component3/component4)
        SymrotangCheck = atan2(component3/component4,component1/component2);
        SymrotCheck = CreateRotationMatrix(0,0,SymrotangCheck,"AZYZ");
        Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what Z rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = Alpha+pi; BetaOut = Beta; GammaOut = Gamma;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZYZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component3 = Arel1(3,2) - Arel1(3,1)*Mcheck(3,2)/Mcheck(3,1);
            component4 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            %SymrotangCheckcos = acos(component1/component2) %a Z rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,0,SymrotangCheck,"AZYZ");
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what Z rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = 2*pi-Alpha; BetaOut = pi-Beta; GammaOut = Gamma+pi;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZYZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component3 = Arel1(3,2) - Arel1(3,1)*Mcheck(3,2)/Mcheck(3,1);
            component4 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            %SymrotangCheckcos = acos(component1/component2) % a Z rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,0,SymrotangCheck,"AZYZ");
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what Z rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = pi-Alpha; BetaOut = pi-Beta; GammaOut = Gamma+pi;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZYZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component3 = Arel1(3,2) - Arel1(3,1)*Mcheck(3,2)/Mcheck(3,1);
            component4 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            %SymrotangCheckcos = acos(component1/component2)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,0,SymrotangCheck,"AZYZ"); % A Z rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what Z rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            disp('Failed isequal check at (tryallanglestest)  Sym1 == 0 && Sym2 == 1, please contact the authors for bughunting')
            AlphaOut = Alpha; BetaOut = Alpha; GammaOut = Alpha;
        end
        
    end
    
    if PASv2(2) == PASv2(3)
        AlphaOut = Alpha; BetaOut = Beta; GammaOut = Gamma;
        RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZYZ");
        Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
        component1 = Arel1(3,1) + Arel1(2,1)*Mcheck(2,1)/Mcheck(3,1);
        component2 = Mcheck(3,1) + Mcheck(2,1)*Mcheck(2,1)/Mcheck(3,1);
        component3 = Arel1(3,1) - Arel1(2,1)*Mcheck(3,1)/Mcheck(2,1);
        component4 = Mcheck(2,1) + Mcheck(3,1)*Mcheck(3,1)/Mcheck(2,1);
        %SymrotangCheckcos = acos(component1/component2)
        %SymrotangChecksin = asin(component3/component4)
        SymrotangCheck = atan2(component3/component4,component1/component2);
        SymrotCheck = CreateRotationMatrix(0,SymrotangCheck,0,"AZXZ"); % X rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
        Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what X rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = Alpha+pi; BetaOut = Beta; GammaOut = Gamma;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZYZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component3 = Arel1(3,1) - Arel1(2,1)*Mcheck(3,1)/Mcheck(2,1);
            component4 = Mcheck(2,1) + Mcheck(3,1)*Mcheck(3,1)/Mcheck(2,1);
            %SymrotangCheckcos = acos(component1/component2)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,SymrotangCheck,0,"AZXZ"); % X rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what X rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = 2*pi-Alpha; BetaOut = pi-Beta; GammaOut = Gamma+pi;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZYZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component3 = Arel1(3,1) - Arel1(2,1)*Mcheck(3,1)/Mcheck(2,1);
            component4 = Mcheck(2,1) + Mcheck(3,1)*Mcheck(3,1)/Mcheck(2,1);
            %SymrotangCheckcos = acos(component1/component2)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,SymrotangCheck,0,"AZXZ"); % X rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what X rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = pi-Alpha; BetaOut = pi-Beta; GammaOut = Gamma+pi;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZYZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component3 = Arel1(3,1) - Arel1(2,1)*Mcheck(3,1)/Mcheck(2,1);
            component4 = Mcheck(2,1) + Mcheck(3,1)*Mcheck(3,1)/Mcheck(2,1);
            %SymrotangCheckcos = acos(component1/component2)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,SymrotangCheck,0,"AZXZ"); % X rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what X rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
            
        else
            disp('Failed isequal check at (tryallanglestest)  Sym1 == 0 && Sym2 == 1, please contact the authors for bughunting')
            AlphaOut = Alpha; BetaOut = Alpha; GammaOut = Alpha;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% ZXZ section

if Mode == "AZXZ" || Mode == "PZXZ"
    if PASv2(1) == PASv2(2)
        AlphaOut = Alpha; BetaOut = Beta; GammaOut = Gamma;
        RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZXZ");
        Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
        component1 = Arel1(3,1) + Arel1(3,2)*Mcheck(3,2)/Mcheck(3,1);
        component2 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
        component3 = Arel1(3,2) - Arel1(3,1)*Mcheck(3,2)/Mcheck(3,1);
        component4 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
        %SymrotangCheckcos = acos(component1/component2) %a Z rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
        %SymrotangChecksin = asin(component3/component4)
        SymrotangCheck = atan2(component3/component4,component1/component2);
        SymrotCheck = CreateRotationMatrix(0,0,SymrotangCheck,"AZYZ");
        Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what Z rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = Alpha+pi; BetaOut = Beta; GammaOut = Gamma;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZXZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component3 = Arel1(3,2) - Arel1(3,1)*Mcheck(3,2)/Mcheck(3,1);
            component4 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            %SymrotangCheckcos = acos(component1/component2) %a Z rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,0,SymrotangCheck,"AZYZ");
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what Z rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = 2*pi-Alpha; BetaOut = pi-Beta; GammaOut = Gamma+pi;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZXZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component3 = Arel1(3,2) - Arel1(3,1)*Mcheck(3,2)/Mcheck(3,1);
            component4 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            %SymrotangCheckcos = acos(component1/component2) % a Z rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,0,SymrotangCheck,"AZYZ");
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what Z rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = pi-Alpha; BetaOut = pi-Beta; GammaOut = Gamma+pi;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZXZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            component3 = Arel1(3,2) - Arel1(3,1)*Mcheck(3,2)/Mcheck(3,1);
            component4 = Mcheck(3,1) + Mcheck(3,2)*Mcheck(3,2)/Mcheck(3,1);
            %SymrotangCheckcos = acos(component1/component2)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,0,SymrotangCheck,"AZYZ"); % A Z rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what Z rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            disp('Failed isequal check at (tryallanglestest)  Sym1 == 0 && Sym2 == 1, please contact the authors for bughunting')
            AlphaOut = Alpha; BetaOut = Alpha; GammaOut = Alpha;
        end
        
    end
    
    if PASv2(2) == PASv2(3)
        AlphaOut = Alpha; BetaOut = Beta; GammaOut = Gamma;
        RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZXZ");
        Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
        component1 = Arel1(3,1) + Arel1(2,1)*Mcheck(2,1)/Mcheck(3,1);
        component2 = Mcheck(3,1) + Mcheck(2,1)*Mcheck(2,1)/Mcheck(3,1);
        component3 = Arel1(3,1) - Arel1(2,1)*Mcheck(3,1)/Mcheck(2,1);
        component4 = Mcheck(2,1) + Mcheck(3,1)*Mcheck(3,1)/Mcheck(2,1);
        %SymrotangCheckcos = acos(component1/component2)
        %SymrotangChecksin = asin(component3/component4)
        SymrotangCheck = atan2(component3/component4,component1/component2);
        SymrotCheck = CreateRotationMatrix(0,SymrotangCheck,0,"AZXZ"); % X rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
        Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what X rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = Alpha+pi; BetaOut = Beta; GammaOut = Gamma;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZXZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component3 = Arel1(3,1) - Arel1(2,1)*Mcheck(3,1)/Mcheck(2,1);
            component4 = Mcheck(2,1) + Mcheck(3,1)*Mcheck(3,1)/Mcheck(2,1);
            %SymrotangCheckcos = acos(component1/component2)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,SymrotangCheck,0,"AZXZ"); % X rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what X rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = 2*pi-Alpha; BetaOut = pi-Beta; GammaOut = Gamma+pi;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZXZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component3 = Arel1(3,1) - Arel1(2,1)*Mcheck(3,1)/Mcheck(2,1);
            component4 = Mcheck(2,1) + Mcheck(3,1)*Mcheck(3,1)/Mcheck(2,1);
            %SymrotangCheckcos = acos(component1/component2)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,SymrotangCheck,0,"AZXZ"); % X rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what X rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
        else
            AlphaOut = pi-Alpha; BetaOut = pi-Beta; GammaOut = Gamma+pi;
            RrelCheck = CreateRotationMatrix(AlphaOut,BetaOut,GammaOut,"AZXZ");
            Mcheck = round(RrelCheck*PAS1*RrelCheck^-1,14);
            component1 = Arel1(3,1) + Arel1(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component2 = Mcheck(3,1) + Mcheck(2,1)*Mcheck(2,1)/Mcheck(3,1);
            component3 = Arel1(3,1) - Arel1(2,1)*Mcheck(3,1)/Mcheck(2,1);
            component4 = Mcheck(2,1) + Mcheck(3,1)*Mcheck(3,1)/Mcheck(2,1);
            %SymrotangCheckcos = acos(component1/component2)
            %SymrotangChecksin = asin(component3/component4)
            SymrotangCheck = atan2(component3/component4,component1/component2);
            SymrotCheck = CreateRotationMatrix(0,SymrotangCheck,0,"AZXZ"); % X rotation around the symmetry axix os the now axially symmetric tensor "tensor 2" (tensor "A" in the paper)
            Mcheck=round(SymrotCheck*Mcheck*SymrotCheck^-1,3); % the rotation checks at what X rotation angle (symmetric rotation around tensor 2) the we can replicate Arel1.
        end
        if isequal(round(Arel1,3),round(Mcheck,3))
            
        else
            disp('Failed isequal check at (tryallanglestest)  Sym1 == 0 && Sym2 == 1, please contact the authors for bughunting')
            AlphaOut = Alpha; BetaOut = Alpha; GammaOut = Alpha;
        end
    end
end