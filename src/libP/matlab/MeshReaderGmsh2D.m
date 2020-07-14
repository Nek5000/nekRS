function [Nv, VX, VY, K, EToV] = MeshReaderGmsh2D(FileName)

  %% function [Nv, VX, VY, K, EToV] = MeshReaderGmsh2D(FileName)
  %% Purpose  : Read in basic grid information to build grid

  file = fopen(FileName, 'rt');
  line = fgetl(file);

  while(line ~= -1)
    
    if strfind(line,'$Nodes')
      
      line = fgetl(file);
      Nv = str2num(line);
      VX = zeros(1,Nv);
      VY = zeros(1,Nv);
      
      for i=1:Nv
	line = fgetl(file);
	out = sscanf(line,'%i%f%f%f');
		     VX(i) = out(2);
		     VY(i) = out(3);
		   end
		   end
		     
		     if strfind(line,'$Elements')
		     
		     line = fgetl(file);
		     K = str2num(line);
		     EToV = zeros(K,3);
		     
		     kTri = 0;
		     for i=1:K
		     line = fgetl(file);
		     out = sscanf(line,'%i%i');
				  if out(2) == 2
				  kTri = kTri+1;
				  out = sscanf(line,'%i%i%i%i%i%i%i%i');
					       EToV(kTri,1) = out(6);
					       EToV(kTri,2) = out(7);
					       EToV(kTri,3) = out(8);
					     end
					     end
					       
					       K = kTri;
					       EToV = EToV(1:K,1:3);
					     end
					       
					       line = fgetl(file);
					     end
