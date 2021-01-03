function DT = CreateDipolar(coorda,coordb)
 Vdip = coordb - coorda;
 r12 = sqrt(sum(Vdip.*Vdip));
 DT = [(3*Vdip(1)^2/(r12^2) - 1), 3*Vdip(1)*Vdip(2)/(r12^2), 3*Vdip(1)*Vdip(3)/(r12^2); 3*Vdip(1)*Vdip(2)/(r12^2),  (3*Vdip(2)^2/(r12^2) - 1), 3*Vdip(2)*Vdip(3)/(r12^2); 3*Vdip(1)*Vdip(3)/(r12^2), 3*Vdip(2)*Vdip(3)/(r12^2), (3*Vdip(3)^2/(r12^2) - 1)];
end