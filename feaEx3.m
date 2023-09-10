%% Material properties and Geometry

E = 200e6;
L= [25*0.0254, 75*0.0254];
w1 = 2*0.0254;
w2 = 5*0.0254;
wavg = (w1+w2)/2;
t = 0.5*0.0254;
A = [w1*t, wavg*t];

%% Meshing

ne = 2;
nn = 3;
nne = 2;
dofn = 1;
dofe = dofn*nne;
tdof = nn*dofn;

Kg = zeros(tdof,tdof);
Fg = zeros(tdof, 1);
q0 = 0;

%% Assembly of Stiffness Matrix

for e=1:ne
    Ke = (E*A(e)/L(e))*[1, -1;-1, 1];
    Fe = (q0*L(e)/2)*[1;1];
    for i=1:dofe
        for j = 1:dofe
            Kg(i+e-1,j+e-1)  = Kg(i+e-1,j+e-1) + Ke(i,j);
        end
        Fg(i+e-1,1) = Fg(i+e-1,1) + Fe(i,1);
    end
end

KgO = Kg;
FgO = Fg;

%% Boundary Conditions

Fg(1,1) = 135;
Kg(:,3) = 0;
Kg(3,:) = 0;
Kg(3,3) = 1;
Fg(3,1) = 0;
Ug = linsolve(Kg, Fg);
Rg = (KgO*Ug)-FgO; 