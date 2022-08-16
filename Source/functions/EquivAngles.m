function [EulerT] = EquivAngles(alpha,beta,gamma,Mode)
%Active rotation solutions: {a,b,g+Pi}, {a+Pi,Pi-b,Pi-g}, and {a+Pi,Pi-b,2Pi-g}
%Passive rotation solutions: {a+Pi,b,g }, {Pi-a,Pi-b,Pi+g}, and {2Pi-a,Pi-b,Pi+g}
if Mode == "AZYZ"
    EulerT = zeros(4,3);
    EulerT(1,:) = [alpha, beta, gamma];
    EulerT(2,:) = [alpha, beta, pi + gamma];
    EulerT(3,:) = [pi + alpha, pi - beta, pi - gamma];
    EulerT(4,:) = [pi + alpha, pi - beta, 2*pi - gamma];
end

if Mode == "PZYZ"
    EulerT = zeros(4,3);
    EulerT(1,:) = [alpha, beta, gamma];
    EulerT(2,:) = [pi + alpha, beta, gamma];
    EulerT(3,:) = [pi - alpha, pi - beta, pi + gamma];
    EulerT(4,:) = [2*pi - alpha, pi - beta,pi + gamma];
end

if Mode == "AZXZ"
    EulerT = zeros(4,3);
    EulerT(1,:) = [alpha, beta, gamma];
    EulerT(2,:) = [alpha, beta, pi + gamma];
    EulerT(3,:) = [pi + alpha, pi - beta, pi - gamma];
    EulerT(4,:) = [pi + alpha, pi - beta, 2*pi - gamma];
end



if Mode == "PZXZ"
    EulerT = zeros(4,3);
    EulerT(1,:) = [alpha, beta, gamma];
    EulerT(2,:) = [pi + alpha, beta, gamma];
    EulerT(3,:) = [pi - alpha, pi - beta, pi + gamma];
    EulerT(4,:) = [2*pi - alpha, pi - beta,pi + gamma];
end

Einda = EulerT(:,1)<0;
Eindb = EulerT(:,2)<0;
Eindg = EulerT(:,3)<0;
if sum(Einda)>0
EulerT(Einda,1) = mod(EulerT(Einda,1),2*pi);
end
if sum(Eindb)>0
EulerT(Eindb,2) = mod(EulerT(Eindb,2),2*pi);
end
if sum(Eindg)>0
EulerT(Eindg,3) = mod(EulerT(Eindg,3),2*pi);
end

Einda = EulerT(:,1)>=2*pi;
Eindb = EulerT(:,2)>=2*pi;
Eindg = EulerT(:,3)>=2*pi;
if sum(Einda)>0
EulerT(Einda,1) = mod(EulerT(Einda,1),2*pi);
end
if sum(Eindb)>0
EulerT(Eindb,2) = mod(EulerT(Eindb,2),2*pi);
end
if sum(Eindg)>0
EulerT(Eindg,3) = mod(EulerT(Eindg,3),2*pi);
end

end