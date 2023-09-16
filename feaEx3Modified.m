E = 200e6;
L = [25*0.0254, 75*0.0254];
ne1 = 2;
ne2 = [2];
nne = 2;
dofn = 1;
dofe = nne*dofn;
t = 0.5*0.0254;
w1 = 2*0.0254;
w2 = 5*0.0254;
wg1 = w1;
Les = 0;
disps = zeros(1, length(ne2));
for e = 1:length(ne2)
    ne = ne1 + ne2(e);
    nn = ne+1;
    tdof = nn*dofn;
    Le = zeros(1, ne);
    Ae = zeros(1, ne);
    for f = 1:ne1
        Le(1,f) = L(1,1)/ne1;
        Ae(1,f) = w1*t;
    end
    for g = (ne1+1):(ne1+ne2(e))
        Le(1,g) = L(1,2)/ne2(e);
        Les = Les + Le(1,g);
        weq = ((3/75)*Les)+(2*0.0254);
        if g==ne1+1
            wg2 = weq;
        else 
            wg1 = wg2;
            wg2 = weq;
        end
        wgAvg = (wg1+wg2)/2;
        Ae(1,g) = wgAvg*t;
    end

    Kg = zeros(tdof, tdof);
    Fg = zeros (tdof, 1);
    
    for o=1:ne
        Ke = (E*Ae(o)/Le(o))*[1, -1;-1, 1];
        q0 = 0;
        Fe = (q0*Le(o)/2)*[1;1];
        for i = 1:dofe  
            for j =1:dofe
                Kg(i+o-1,j+o-1) = Kg(i+o-1,j+o-1) + Ke(i,j);
            end
            Fg(i+o-1, 1) = Fg(i+o-1, 1) + Fe(i,1);
        end
    end
    
    Fg(1,1) = -135;
    Kg(:,tdof) = 0;
    Kg(tdof,:) = 0;
    Kg(tdof,tdof) = 1;
    Fg(tdof,1) = 0;
    Ug = linsolve(Kg, Fg);

    plot(ne,Ug(1,1),'r*--');
    hold on;
    disps(1,e) = Ug(1,1);

    Les = 0;
end