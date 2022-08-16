function [EulerT] = EquivRelAngles(alpha,beta,gamma,Mode)
if Mode == "AZYZ"
    EulerT = zeros(16,3);
    EulerT(1,:) = [alpha, beta, gamma];
    EulerT(2,:) = [alpha + pi, pi - beta, 2*pi - gamma];
    EulerT(3,:) = [alpha + pi, pi - beta, pi - gamma];
    EulerT(4,:) = [alpha, beta, pi + gamma];
    
    EulerT(5,:) = [2*pi - alpha, pi - beta, pi + gamma];
    EulerT(6,:) = [pi - alpha, beta, pi - gamma];
    EulerT(7,:) = [pi - alpha, beta, 2*pi - gamma];
    EulerT(8,:) = [2*pi - alpha, pi - beta, gamma];
    
    EulerT(9,:) = [pi - alpha, pi - beta, pi + gamma];
    EulerT(10,:) = [2*pi - alpha, beta, pi - gamma];
    EulerT(11,:) = [2*pi - alpha, beta,2*pi - gamma];
    EulerT(12,:) = [pi - alpha, pi - beta, gamma];
    
    EulerT(13,:) = [alpha + pi, beta, gamma];
    EulerT(14,:) = [alpha, pi - beta, 2*pi - gamma];
    EulerT(15,:) = [alpha, pi - beta, pi - gamma];
    EulerT(16,:) = [alpha + pi, beta, pi + gamma];
end

if Mode == "PZYZ"
    EulerT = zeros(16,3);
    EulerT(1,:) = [alpha, beta, gamma];
    EulerT(2,:) = [2*pi - alpha, pi - beta, pi + gamma];
    EulerT(3,:) = [pi - alpha, pi - beta, pi + gamma];
    EulerT(4,:) = [pi + alpha, beta, gamma];
    
    EulerT(5,:) = [pi + alpha, pi - beta, 2*pi - gamma];
    EulerT(6,:) = [pi - alpha, beta, pi - gamma];
    EulerT(7,:) = [2*pi - alpha, beta, pi - gamma];
    EulerT(8,:) = [alpha, pi - beta, 2*pi - gamma];
    
    EulerT(9,:) = [pi + alpha, pi - beta, pi - gamma];
    EulerT(10,:) = [pi - alpha, beta, 2*pi - gamma];
    EulerT(11,:) = [2*pi - alpha, beta,2*pi - gamma];
    EulerT(12,:) = [alpha, pi - beta, pi - gamma];
    
    EulerT(13,:) = [alpha, beta, pi + gamma];
    EulerT(14,:) = [2*pi - alpha, pi - beta, gamma];
    EulerT(15,:) = [pi - alpha, pi - beta, gamma];
    EulerT(16,:) = [pi + alpha, beta, pi + gamma];
end

if Mode == "AZXZ"
    EulerT = zeros(16,3);
    EulerT(1,:) = [alpha, beta, gamma];
    EulerT(2,:) = [pi + alpha, pi - beta,2*pi - gamma];
    EulerT(3,:) = [pi + alpha, pi - beta, pi - gamma];
    EulerT(4,:) = [alpha, beta, pi + gamma];
    
    EulerT(5,:) = [2*pi - alpha,pi - beta, pi + gamma];
    EulerT(6,:) = [pi - alpha, beta, pi - gamma];
    EulerT(7,:) = [pi - alpha, beta,2*pi - gamma];
    EulerT(8,:) = [2*pi - alpha,pi - beta, gamma];
    
    EulerT(9,:) = [pi - alpha,pi - beta, pi + gamma];
    EulerT(10,:) = [2*pi - alpha, beta, pi - gamma];
    EulerT(11,:) = [2*pi - alpha, beta,2*pi - gamma];
    EulerT(12,:) = [pi - alpha,pi - beta, gamma];
    
    EulerT(13,:) = [pi + alpha, beta, gamma];
    EulerT(14,:) = [alpha,pi - beta,2*pi - gamma];
    EulerT(15,:) = [alpha,pi - beta,pi - gamma];
    EulerT(16,:) = [pi + alpha, beta,pi + gamma];
end



if Mode == "PZXZ"
    EulerT = zeros(16,3);
    EulerT(1,:) = [alpha,beta, gamma];
    EulerT(2,:) = [2*pi - alpha,pi - beta, pi + gamma];
    EulerT(3,:) = [pi - alpha,pi - beta, pi + gamma];
    EulerT(4,:) = [pi + alpha, beta, gamma];
    
    EulerT(5,:) = [pi + alpha,pi - beta,2*pi - gamma];
    EulerT(6,:) = [pi - alpha, beta,pi - gamma];
    EulerT(7,:) = [2*pi - alpha, beta,pi - gamma];
    EulerT(8,:) = [alpha,pi - beta,2*pi - gamma];
    
    EulerT(9,:) = [pi + alpha,pi - beta,pi - gamma];
    EulerT(10,:) = [pi - alpha, beta,2*pi - gamma];
    EulerT(11,:) = [2*pi - alpha, beta,2*pi - gamma];
    EulerT(12,:) = [alpha,pi - beta,pi - gamma];
    
    EulerT(13,:) = [alpha, beta, pi + gamma];
    EulerT(14,:) = [2*pi - alpha,pi - beta, gamma];
    EulerT(15,:) = [pi - alpha, pi - beta, gamma];
    EulerT(16,:) = [pi + alpha, beta,pi + gamma];
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