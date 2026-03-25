nBarras=60;
numNodes=25;
conectividades=[1     2
     1     3
     1     4
     1     5
     1     6
     1     7
     2     3
     3     4
     4     5
     5     6
     6     7
     7     2
     2     8
     2     9
     3     9
     3    10
     4    10
     4    11
     5    11
     5    12
     6    12
     6    13
     7    13
     7     8
     8     9
     9    10
    10    11
    11    12
    12    13
    13     8
     8    25
     8    15
     9    15
     9    17
    10    17
    10    19
    11    19
    11    21
    12    21
    12    23
    13    23
    13    25
    15    17
    17    19
    19    21
    21    23
    23    25
    25    15
    15    14
    15    16
    17    16
    17    18
    19    18
    19    20
    21    20
    21    22
    23    22
    23    24
    25    24
    25    14];           
numColConect=size(conectividades,2);
conectVetor=[];
for cLncv=1:size(conectividades,1)
    conectVetor=[conectVetor conectividades(cLncv,:)];
end    
nodesC = zeros(1,numNodes*3);
nodesC = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0]; 

numConst=size(find(nodesC==1),2);

MatAux1=[1       2       -10000]; 
     
LinMatAux1=size(MatAux1,1);
nodesL=zeros(1,3*numNodes);
for cLinMatAux1=1:LinMatAux1
    nodesL(MatAux1(cLinMatAux1,1)*3-(3-(MatAux1(cLinMatAux1,2))))=MatAux1(cLinMatAux1,3);
end
CoordAux=[0   34.7070         0
  -21.6506   32.7070  -12.5000
         0   32.7070  -25.0000
   21.6506   32.7070  -12.5000
   21.6506   32.7070   12.5000
         0   32.7070   25.0000
  -21.6506   32.7070   12.5000
  -50.0000   26.4910         0
  -25.0000   26.4910  -43.3013
   25.0000   26.4910  -43.3013
   50.0000   26.4910         0
   25.0000   26.4910   43.3013
  -25.0000   26.4910   43.3013
 -100.0000         0         0
  -64.9519   14.1640  -37.5000
  -50.0000         0  -86.6025
         0   14.1640  -75.0000
   50.0000         0  -86.6025
   64.9519   14.1640  -37.5000
  100.0000         0         0
   64.9519   14.1640   37.5000
   50.0000         0   86.6025
         0   14.1640   75.0000
  -50.0000         0   86.6025
  -64.9519   14.1640   37.5000]./10;
CoordVetor=[];
for cLcv=1:size(CoordAux,1)
    CoordVetor=[CoordVetor CoordAux(cLcv,:)];
end       
for aaaaux=1:nBarras
    elasticity(1,aaaaux)=73000000000;
end
ro = 2770;

nodesCC=[14 16 18 20 22 24];
tamNodesCC=size(nodesCC,2);

aalfa=zeros(1,nBarras);

