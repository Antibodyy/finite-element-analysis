E = 200e9;
A = 3250e-6;
L1 = 3.6;
L2 = 3.118;

ne = 11;
nn = 7;
nne = 2;
dofn = 2;
dofe = nne*dofn;
tdof = nn*dofn;

nodalCONN = [1 2;1 3;2 3;2 4;3 4;3 5;4 5;4 6;5 6;5 7;6 7];
CONN = [1 2 3 4;1 2 5 6;3 4 5 6;3 4 7 8;5 6 7 8;5 6 9 10;7 8 9 10;7 8 11 12;9 10 11 12;9 10 13 14;11 12 13 14];

x = [0 L1/2 L1 3*L1/2 2*L1 5*L1/2 3*L1];
y = [0 L2 0 L2 0 L2 0];

for i = 1:ne
    Le(i) = sqrt((y(nodalCONN(i,2))-y(nodalCONN(i,1)))^2+(x(nodalCONN(i,2))-x(nodalCONN(i,1)))^2);
    theta(i) = atan2d((y(nodalCONN(i,2))-y(nodalCONN(i,1))),(x(nodalCONN(i,2))-x(nodalCONN(i,1))));
end

%% Assembly

Kg = zeros(tdof,tdof);
Fg = zeros(tdof,1);

for e=1:ne
    KeL = (E*A/Le(e))*[1 -1;-1 1];
    Te = [cos(theta(e)*pi/180) sin(theta(e)*pi/180) 0 0;0 0 cos(theta(e)*pi/180) sin(theta(e)*pi/180)];
    KeG = transpose(Te)*KeL*Te;
    for i = 1:dofe
        for j = 1:dofe
            Kg(CONN(e,i),CONN(e,j)) = Kg(CONN(e,i),CONN(e,j)) + KeG(i,j);
        end
    end
end

%% Boundary Conditions

Kg(1,:) = 0;
Kg(:,1) = 0;
Kg(2,:) = 0;
Kg(:,2) = 0;
Kg(14,:) = 0;
Fg(2,1) = -280e3;
Fg(6,1) = -210e3;
Fg(10,1) = -280e3;
Fg(14,1) = -360e3;
Fg(2,1) = 0;
Fg(1,1) = 0;
Fg(14,1) = 0;
Kg(1,1) = 1;
Kg(2,2) = 1;
Kg(14,14) = 1;

Ug = linsolve(Kg,Fg);