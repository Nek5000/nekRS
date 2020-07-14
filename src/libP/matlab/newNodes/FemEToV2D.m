function [refEToV] = FemEToV2D(N,req,seq,nodeType)

refEToV = delaunayOriented2D(req',seq');
if strcmp(nodeType,'EI')
    if (N==2)
        %refEToV(2,:) = [4,1,2];
        %refEToV(5,:) = [4,2,7];
    end
end

if strcmp(nodeType,'SW')
    %{
    if (N==6)
        refEToV(2,:) = [29,19,30];
        refEToV(40,:) = [30,28,29];
    elseif (N==4)
        refEToV(6,:) = [13,18,17];
        refEToV(13,:) = [15,19,17];
        refEToV(17,:) = [17,18,16];
        refEToV(10,:) = [19,13,17];
    elseif (N==8)
        refEToV(56,:) = [25,41,43];
        refEToV(60,:) = [25,42,41];
        refEToV(19,:) = [40,41,42];
        refEToV(58,:) = [43,41,39];
    end 
    %}
end

end