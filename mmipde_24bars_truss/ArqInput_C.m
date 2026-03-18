nBarras=24;
numNodes=13;
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
                3 8;
                4 8;
                6 9;
                7 9;
                2 12;
                3 12;
                2 13;
                7 13;
                4 10;
                5 10;
                5 11;
                6 11];           
numColConect=size(conectividades,2);
conectVetor=[];
for cLncv=1:size(conectividades,1)
    conectVetor=[conectVetor conectividades(cLncv,:)];
end    
nodesC = zeros(1,numNodes*3);
nodesC = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; 

numConst=size(find(nodesC==1),2);

MatAux1=[1      2       -50000]; 
     
LinMatAux1=size(MatAux1,1);
nodesL=zeros(1,3*numNodes);
for cLinMatAux1=1:LinMatAux1
    nodesL(MatAux1(cLinMatAux1,1)*3-(3-(MatAux1(cLinMatAux1,2))))=MatAux1(cLinMatAux1,3);
end
CoordAux=[0    3.0000         0
   -9.0000    1.5000         0
   -4.5000    1.5000    7.7942
    4.5000    1.5000    7.7942
    9.0000    1.5000         0
    4.5000    1.5000   -7.7942
   -4.5000    1.5000   -7.7942
         0         0   15.5885
         0         0  -15.5885
   13.5000         0    7.7942
   13.5000         0   -7.7942
  -13.5000         0    7.7942
  -13.5000         0   -7.7942];
CoordVetor=[];
for cLcv=1:size(CoordAux,1)
    CoordVetor=[CoordVetor CoordAux(cLcv,:)];
end       
for aaaaux=1:nBarras
    elasticity(1,aaaaux)=73000000000;
end
ro = 2770;

nodesCC=[8:13];
tamNodesCC=size(nodesCC,2);

aalfa=zeros(1,nBarras);

