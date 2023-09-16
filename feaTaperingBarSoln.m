
E=200e6; %Young's Modulus
W1=5*.0254; %Larger width
W2=2*0.0254; %Smaller width
t=0.5*0.0254; %Thickness
L1=75*0.0254; %Length of tapered section
L2=25*0.0254; %Length of uniform section
A = [(W1+W2)*t/2 W2*t];
F=-135e3; %Concentrated load


Lt=1.905; %Length of tapered section
L2 = L2/2;
for ne=[4 8 16 32]
W1=5*0.0254;
W2=2*0.0254;
nne=2; %Number of nodes per element
nn=ne*nne-(ne-1); %Number of nodes
dofn=1; %Degrees of freedom per node
dofe=dofn*nne; %Degrees of freedom per element
tdof=dofn*nn; %Total Degrees of freedom
L1=Lt/(ne-2); %Length of tapered element element
slope=(W1-W2)/Lt; %Change in width per unit length
CONN= zeros(ne,2); %Connectivity matrix
KG=zeros(tdof,tdof); %Global stiffness matrix
FGC=zeros(tdof,1); %Global load vector for concentrated load
for i=1:ne
CONN(i,:)=[i,i+1]; %Nodal Connectivity matrix
W2=W1-slope*L1;
Wavg=(W1+W2)/2; %Average width
if i<=ne-2
A=Wavg*t; %Area of tapered cross-section
Le = L1;
else
A = 0.0508*t;
Le = L2;
end
Ke=E*A/Le*[1 -1; -1 1]; %Elemental stiffness matrix
W1=W2;
for j=1:dofe
for k=1:dofe
KG(CONN(i,j),CONN(i,k))= KG(CONN(i,j),CONN(i,k))+Ke(j,k);
end
end
end
FGC(tdof,1)=F;
FG = FGC;
KGR=KG;
FGR=FG;
KG(1,:)=[];
KG(:,1)=[];
FG(1,:)=[];
UG= linsolve(KG,FG);
Disp_at_A=UG(tdof-1,1);
plot(ne,Disp_at_A,'r*--')
hold on
end
xticks([4 8 16 32]);
xlabel('Number of elements');
ylabel('Deflection at the free end(in Metre)');