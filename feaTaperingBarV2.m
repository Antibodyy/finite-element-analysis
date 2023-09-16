E = 200e6;
L = [25*0.0254 75*0.0254];
t = 0.5*0.0254;
w1 = 2*0.0254;
wg1 = w1;

ne1 = 2;
ne2 = [2, 6, 14, 32];

for e = 1:length(ne2)
    LeS = 0;
    ne = ne1 + ne2(e);
    nn = ne+1;
    nne = 2;
    dofn = 1;
    tdof = nn*dofn;
    dofe = nne*dofn;
    Le = zeros(ne,1);
    Ae = zeros(ne,1);
    
    for a = 1:ne1
        Le(a,1) = L(1,1)/ne1;
        Ae(a,1) = t*w1;
    end

    for b = ne1+1:ne
        Le(b,1) = L(1,2)/ne2(e);
        LeS = Le(b,1) + LeS;
        wg2 = ((3/75)*LeS)+(2*0.0254);
        wAvg = (wg1+wg2)/2;
        Ae(b,1) = wAvg*t;
        wg1 = wg2;
    end

    Kg = zeros(tdof,tdof);
    Fg = zeros(tdof,1);
    q0 = 0;

    for i = 1:ne
        Ke = (E*Ae(i,1)/Le(i,1))*[1 -1;-1 1];
        Fe = (q0*Le(i,1)/2)*[1;1];
        for j = 1:dofe
            for k = 1:dofe
                Kg(i+j-1,i+k-1) = Kg(i+j-1,i+k-1) + Ke(j,k);
            end
            Fg(i+j-1,1) = Fg(i+j-1,1) + Fe(j,1);

        end
    end

    Fg(1,1) = -135;
    Kg(nn,:) = 0;
    Kg(:,nn) = 0;
    Kg(nn,nn) = 1;

    Ug = linsolve(Kg,Fg);
    plot(ne, Ug(1,1),'r*--');
    hold on;

end


