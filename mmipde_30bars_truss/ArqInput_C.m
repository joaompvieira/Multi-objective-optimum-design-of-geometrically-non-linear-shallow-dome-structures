nBarras=30;
numNodes=19;
conectividades=[1 2;
                1 3;
                1 4;
                1 5;
                1 6;
                1 7;
                2 3;
                3 4;
                4 5;
                5 6;
                6 7;
                7 2;
                2 9
                3 11
                4 13
                5 15
                6 17
                7 19
                2 8
                2 10
                3 10
                3 12
                4 12
                4 14
                5 14
                5 16
                6 16
                6 18
                7 18
                7 8];           
numColConect=size(conectividades,2);
conectVetor=[];
for cLncv=1:size(conectividades,1)
    conectVetor=[conectVetor conectividades(cLncv,:)];
end    
nodesC = zeros(1,numNodes*3);
nodesC = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; 

numConst=size(find(nodesC==1),2);

MatAux1=[1      2       -2000]; 
     
LinMatAux1=size(MatAux1,1);
nodesL=zeros(1,3*numNodes);
for cLinMatAux1=1:LinMatAux1
    nodesL(MatAux1(cLinMatAux1,1)*3-(3-(MatAux1(cLinMatAux1,2))))=MatAux1(cLinMatAux1,3);
end
CoordAux=[0.000  0.000  85.912
0.000  -360  64.662
311.769  -180  64.662
311.769  180  64.662
0.000  360  64.662
-311.769  180  64.662
-311.769  -180  64.662
-311.769  -540  21.709
0.000 -720  0.000
311.769  -540  21.709
623.538 -360  0.000
623.538  0.000  21.709
623.538 360  0.000
311.769  540  21.709
0.000 720  0.000
-311.769  540  21.709
-623.538 360  0.000
-623.538  0.000  21.709
-623.538 -360  0.000];
Coordn = CoordAux(:,2);
CoordAux(:,2) = CoordAux(:,3);
CoordAux(:,3) = Coordn(:,1);
CoordVetor=[];
for cLcv=1:size(CoordAux,1)
    CoordVetor=[CoordVetor CoordAux(cLcv,:)];
end       
for aaaaux=1:nBarras
    elasticity(1,aaaaux)=10000000;
end
ro = 0.1;

nodesCC=[8:19];
tamNodesCC=size(nodesCC,2);

aalfa=zeros(1,nBarras);

