function [r,s] = NewEquiNodes2D(N,nodeType)

if (strcmp(nodeType,'WB'))
    Nbndry = N;
    Ninter = N;
    Layered = 0;
elseif (strcmp(nodeType,'EI')) 
    Nbndry = N;
    Ninter = N+1;
    Layered = 0;
elseif (strcmp(nodeType,'SW')) 
    Nbndry = N;
    Ninter = N+1;
    Layered = 1;
else
    error('bad node type in NewEuiNodes2D')
end

if Layered
    step = Ninter-Nbndry;
    M = Nbndry;
    r = [];
    s = [];
    scale = 1;
    d = .1/(Ninter*Ninter);
    while (M>-1) 
        if M==0
            rM = 0;
            sM = 0;
            ids = 1;
        else
            [rM,sM] = EquiNodes2D(M);
            %extract the bounary nodes
            ids = find((sM+1/sqrt(3))-(rM+1)*3/sqrt(3) > -d | sM+1/sqrt(3) < d | (sM-2/sqrt(3))+rM*3/sqrt(3) > -d);
        end

        r = [r;scale*rM(ids)];
        s = [s;scale*sM(ids)];

        scale = scale*(1-3/(M+step));

        M = M-3+step;
    end
else
    d = .1/(Ninter*Ninter);
    [r1,s1] = EquiNodes2D(Ninter);
    
    ids = find((s1+1/sqrt(3))-(r1+1)*3/sqrt(3) < -d & s1+1/sqrt(3) > d & (s1-2/sqrt(3))+r1*3/sqrt(3) < -d);
    r1 = r1(ids);
    s1 = s1(ids);

    [r,s] = EquiNodes2D(Nbndry);
    ids = find((s+1/sqrt(3))-(r+1)*3/sqrt(3) > -d | s+1/sqrt(3) < d | (s-2/sqrt(3))+r*3/sqrt(3) > -d);
    r = [r(ids);r1];
    s = [s(ids);s1];
end

end