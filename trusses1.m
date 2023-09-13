%% Material Properties and Meshing

nn = 3;
ne = 3;
tdof = nn*2;
dofe = 4;
dofn = 2;
E = 200e9;
A = 30e-4;

nodalCONN = [1 2;2 3;1 3];
CONN = [1 2 3 4;3 4 5 6;1 2 5 6];
x = [0 100e-2 100e-2];
y = [0 0 100e-2];
for i = 1:ne
    L(i) = sqrt(((y(nodalCONN(i,2))-y(nodalCONN(i,1)))^2)+(x(nodalCONN(i,2))-x(nodalCONN(i,1)))^2);
    theta(i) = atan2d((y(nodalCONN(i,2))-y(nodalCONN(i,1))), (x(nodalCONN(i,2))-x(nodalCONN(i,1))));
end

%% Assembly

KG = zeros(tdof,tdof);
FG = zeros(tdof,1);

for e = 1:ne
    KeL = (E*A/L(e))*[1 -1;-1 1];
    Te = [cos(theta(e)*pi/180) sin(theta(e)*pi/180) 0 0;0 0 cos(theta(e)*pi/180) sin(theta(e)*pi/180)];
    KeG = transpose(Te)*KeL*Te;
    for i = 1:dofe
        for j = 1:dofe
            KG(CONN(e,i),CONN(e,j)) = KG(CONN(e,i),CONN(e,j)) + KeG(i,j);
        end
    end
end

%% Boundary Conditions

KG(1,:) = 0;
KG(:,1) = 0;
KG(2,:) = 0;
KG(:,2) = 0;
KG(3,:) = 0;
KG(:,3) = 0;
KG(4,:) = 0;
KG(:,4) = 0;
FG(5,1) = 100e3;
FG(6,1) = -200e3;
KG(1,1) = 1;
KG(2,2) = 1;
KG(3,3) = 1;
KG(4,4) = 1;

UG = linsolve(KG,FG);