%% Material Properties and Geometry Definition

E = 200e9;
I = 4000000e-12;
L = [1 1];
q0 = [0 -12e3];

%% Meshing

nn = 3;
ne = 2;
nne = 2;
dofn = 2;
dofe = dofn*nne;
tdof = nn*dofn;

Kg = zeros(tdof,tdof);
Fg = zeros(tdof,1);

CONN = [1 2 3 4;3 4 5 6];

%% Assembly of Stiffness Matrix

for e = 1:ne  
    Ke = (E*I/(L(e)^3))*[12 6*L(e) -12 6*L(e); 6*L(e) 4*(L(e)^2) -6*L(e) 2*(L(e)^2); -12 -6*L(e) 12 -6*L(e); 6*L(e) 2*(L(e)^2) -6*L(e) 4*(L(e)^2)];
    Fe = [q0(e)*L(e)/2; q0(e)*(L(e)^2)/12; q0(e)*L(e)/2; -q0(e)*(L(e)^2)/12];
    for i = 1:dofe
        for j = 1:dofe
            Kg(CONN(e,i),CONN(e,j)) = Kg(CONN(e,i),CONN(e,j)) + Ke(i,j);
        end
        Fg(CONN(e,i),1) = Fg(CONN(e,i),1) + Fe(i,1);
    end
end

KgO = Kg;

%% Boundary Conditions

Kg(1,:) = 0;
Kg(:,1) = 0;
Kg(1,1) = 1;
Fg(1,1) = 0;
Kg(2,:) = 0;
Kg(:,2) = 0;
Kg(2,2) = 1;
Fg(2,1) = 0;
Kg(3,:) = 0;
Kg(:,3) = 0;
Kg(3,3) = 1;
Fg(3,1) = 0;
Kg(5,:) = 0;
Kg(:,5) = 0;
Kg(5,5) = 1;
Fg(5,1) = 0;

wg = linsolve(Kg,Fg);
FgO = KgO*wg;
