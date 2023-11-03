%% Geometric and Material Properties
E = 200e9;
mu = 0.3;
Le = 2;
He = 2;
t = 0.02;

%% Mesh Details
nn = 4;
ne = 2;
nne = 3;
dofn = 2;
dofe = nne*dofn;
tdof = nn*dofn;

% Nodal conn. and corresponding connectivity matrix generation 
nodalCONN = [4 1 2;4 2 3];
CONN = zeros(ne,2*nne);
for v = 1:ne
    CONN(v,1) = (nodalCONN(v,1)*2)-1;
    CONN(v,2) = (nodalCONN(v,1)*2);
    CONN(v,3) = (nodalCONN(v,2)*2)-1;
    CONN(v,4) = (nodalCONN(v,2)*2);
    CONN(v,5) = (nodalCONN(v,3)*2)-1;
    CONN(v,6) = (nodalCONN(v,3)*2);
end

coords = zeros(4,2);
coords(1,1) = 2;
coords(1,2) = 0;
coords(2,1) = 2;
coords(2,2) = 2;
coords(3,1) = 0;
coords(3,2) = 2;
coords(4,1) = 0;
coords(4,2) = 0;

%% CST Element Stiffness Matrix Calculation

Kg = zeros(tdof,tdof);
Fg = zeros(tdof,1);

for x = 1:ne
% Elementwise local node definition
    x1 = coords(nodalCONN(x,1),1);
    y1 = coords(nodalCONN(x,1),2);
    x2 = coords(nodalCONN(x,2),1);
    y2 = coords(nodalCONN(x,2),2);
    x3 = coords(nodalCONN(x,3),1);
    y3 = coords(nodalCONN(x,3),2);

    J = [(x1-x3) (x2-x3);(y1-y3) (y2-y3)]; %Jacobian
    A = abs(det(J))/2; %Area of Element
    D = (E/(1-mu^2))*[1 mu 0;mu 1 0;0 0 (1-mu)/2]; %Constant Matrix for Homogenous Material 
    Be = (1/abs(det(J)))*[(y2-y3) 0 (y3-y1) 0 (y1-y2) 0;
                      0 (x3-x2) 0 (x1-x3) 0 (x2-x1);
                      (x3-x2) (y2-y3) (x1-x3) (y3-y1) (x2-x1) (y1-y2)]; 
    Ke = (transpose(Be)*D*Be)*A*t;
    
    for i = 1:dofe
        for j = 1:dofe
            Kg(CONN(x,i),CONN(x,j)) = Kg(CONN(x,i),CONN(x,j)) + Ke(i,j);
        end
    end

    
Fg(4,1) = -30e3;
KgO = Kg;

%% Boundary Conditions

Kg(7,:) = 0;
Kg(:,7)= 0;
Kg(7,7) = 1;
Kg(8,:) = 0;
Kg(:,8)= 0;
Kg(8,8) = 1;
Kg(5,:) = 0;
Kg(:,5)= 0;
Kg(5,5) = 1;
Kg(6,:) = 0;
Kg(:,6)= 0;
Kg(6,6) = 1;
Fg(7,1) = 0;
Fg(8,1) = 0;
Fg(5,1) = 0;
Fg(6,1) = 0;

%% Deflection

ug = linsolve(Kg,Fg);
Rg = KgO*ug;
end