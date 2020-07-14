function [r,s] = NewNodes2D(N,nodeType)

alphaWBLebN = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 ...
                1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258];
alphaWBKappaN = [5.0000 5.0000 3.0902 1.4980 0.0000 0.8985 1.1027 1.4207 1.5333 1.6003 1.5502 1.6045];
alphaWBLebNp1 = [3.4069 3.9919 1.2717 2.2699 1.9084 1.7663 1.6578 1.5990 1.6123 1.6111 1.7392 1.9327];
alphaWBKappaNp1 = [5.0000 5.0000 0.0000 3.2087 2.7064 2.6906 2.7480 2.6873 2.6226 2.5785 2.5150 1.6471];
alphaEILebN = [5.0000 1.9543 0.0005 1.2966 2.0677 1.4920 1.5477 1.4305 1.5795 1.5529 1.7585 1.9543];
alphaEIKappaN = [5.0000 5.0000 4.7045 4.0979 3.2906 3.0382 2.7225 2.6790 2.6684 2.7694 2.8031 2.7811];
alphaEILebNp1 = [1.8280 2.9703 0.0352 0.9370 1.3929 0.7571 1.2899 1.2518 2.0003 1.8798 1.5833 1.8785];
alphaEIKappaNp1 = [5.0000 5.0000 0.0001 2.2073 2.5259 2.7113 2.4368 2.4564 2.3948 2.4346 2.4653 2.4691];
alphaSWLebN = [5.0000 2.5543 0.0006 0.4892 2.5395 2.0847 1.9457 1.8606 1.9219 2.2728 2.1161 2.1146];
alphaSWKappaN = [5.0000 5.0000 4.7045 4.1389 3.5038 3.0405 2.7641 2.7182 2.6271 2.6612 2.6271 2.5668];
alphaSWLebNp1 = [2.9445 3.9913 0.0636 0.0684 2.1170 2.1709 2.1688 2.2124 2.0804 2.0030 2.1696 1.9768];
alphaSWKappaNp1 = [5.0000 5.0000 0.0001 2.5256 2.5295 2.5716 2.5202 2.4872 2.4532 2.4188 2.4005 2.3759];                

if (strcmp(nodeType,'WB')||strcmp(nodeType,'WBLebN'))
    [r,s] = NewEquiNodes2D(N,'WB');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N, r, s, alphaWBLebN(N));
elseif (strcmp(nodeType,'WBKappaN')) 
    [r,s] = NewEquiNodes2D(N,'WB');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N, r, s, alphaWBKappaN(N));
elseif (strcmp(nodeType,'WBLebNp1'))
    [r,s] = NewEquiNodes2D(N+1,'WB');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N+1, r, s, alphaWBLebNp1(N));
elseif (strcmp(nodeType,'WBKappaNp1')) 
    [r,s] = NewEquiNodes2D(N+1,'WB');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N+1, r, s, alphaWBKappaNp1(N));
elseif (strcmp(nodeType,'EILebN')) 
    [r,s] = NewEquiNodes2D(N,'EI');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N, r, s, alphaEILebN(N));
elseif (strcmp(nodeType,'EILebNp1'))
    [r,s] = NewEquiNodes2D(N+1,'EI');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N+1, r, s, alphaEILebNp1(N));
elseif (strcmp(nodeType,'EIKappaN')) 
    [r,s] = NewEquiNodes2D(N,'EI');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N, r, s, alphaEIKappaN(N));
elseif (strcmp(nodeType,'EIKappaNp1'))
    [r,s] = NewEquiNodes2D(N+1,'EI');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N+1, r, s, alphaEIKappaNp1(N));
elseif (strcmp(nodeType,'SWLebN')) 
    [r,s] = NewEquiNodes2D(N,'SW');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N, r, s, alphaSWLebN(N));
elseif (strcmp(nodeType,'SWLebNp1'))
    [r,s] = NewEquiNodes2D(N+1,'SW');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N+1, r, s, alphaSWLebNp1(N));
elseif (strcmp(nodeType,'SWKappaN')) 
    [r,s] = NewEquiNodes2D(N,'SW');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N, r, s, alphaSWKappaN(N));
elseif (strcmp(nodeType,'SWKappaNp1'))
    [r,s] = NewEquiNodes2D(N+1,'SW');
    [r,s] = xytors(r,s);
    [r,s] = WarpBlendTransform2D(N+1, r, s, alphaSWKappaNp1(N));
end

end