%% Material Properties & Geometry
E = [70e9, 200e9];
A = [2400e-6 600e-6];
L = [0.3 0.4];

%% Meshing
ne = 2;
nne = 2;
nn = ne + 1; %May change
dofn = 1;
dofe = nne*dofn;
tdof = dofn*nn;
Kg = zeros(tdof,tdof);
Fg = zeros(tdof,1);
q0 = 0;

%% Assembly of Stiffness Matrix
for e=1:ne
    Ke = ((E(e)*A(e))/L(e))*[1, -1;-1, 1];
    Fe = ((q0*L(e))/2)*[1;1];
    for i = 1:dofe
        for j = 1:dofe
            Kg(i+e-1,j+e-1) = Kg(i+e-1,j+e-1)+Ke(i,j);
        end
        Fg(i+e-1,1) = Fg(i+e-1,1) + Fe(i,1);
    end
end
KgO = Kg;
FgO = Fg;

%% Boundary Conditions

Fg(2,1) = 200e3;
Kg(:,1) = 0;
Kg(:,3) = 0;
Kg(1,:) = 0;
Kg(3,:) = 0;
Kg(1,1) = 1;
Kg(3,3) = 1;
Fg(1,1) = 0;
Fg(3,1) = 0;
Ug = linsolve(Kg,Fg);