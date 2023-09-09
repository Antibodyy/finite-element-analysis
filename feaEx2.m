clc
clear

%% Material Properties and Geometry 
E = [100e9, 70e9];
A = [314.16e-6, 78.54e-6];
L = [1.8, 0.8];

%% Meshing
ne = 2;
nn = ne+1;
nne = 2;
dofn = 1;
dofe = dofn*nne;
tdof = nn*dofn;

Kg = zeros(tdof,tdof);
Fg = zeros(tdof,1);
q0 = 0;

%% Assembly of Stiffness Matrix

for e=1:ne
    Ke = ((E(e)*A(e))/L(e))*[1, -1;-1, 1];
    Fe = (q0*L(e)/2)*[1;1];
    for i = 1:dofe
        for j = 1:dofe
            Kg(i+e-1,j+e-1) = Kg(i+e-1,j+e-1) + Ke(i,j);
        end
        Fg(i+e-1,1) = Fg(i+e-1,1) + Fe(i,1);
    end
end

FgO = Fg;
KgO = Kg;

%% Boundary Conditions

Fg(3,1) = 50e3;
Kg(1,:) = 0;
Kg(:,1) = 0;
Kg(1,1) = 1;
Fg(1,1) = 0;

Ug = linsolve(Kg, Fg);
Exx = (1/L(1,1))*[-1,1]*[Ug(1,1); Ug(2,1)];
Stress900mm = E(1,1)*Exx;
Rg = KgO*Ug-Fg;

