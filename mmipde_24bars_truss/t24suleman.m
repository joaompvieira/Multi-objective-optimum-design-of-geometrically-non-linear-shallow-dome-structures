function [f,g] = t24suleman(Solucao_candidata)    

    Solucao_candidata=Solucao_candidata';
    dados_do_problema;

%     Solucao_candidata = [0.00026544 0.00017481 0.00009082];
    areaSolucao_candidata(1)=Solucao_candidata(1,1);        
    areaSolucao_candidata(2)=Solucao_candidata(1,1);
    areaSolucao_candidata(3)=Solucao_candidata(1,1);
    areaSolucao_candidata(4)=Solucao_candidata(1,1);
    areaSolucao_candidata(5)=Solucao_candidata(1,1);
    areaSolucao_candidata(6)=Solucao_candidata(1,1);
    areaSolucao_candidata(7)=Solucao_candidata(1,2);
    areaSolucao_candidata(8)=Solucao_candidata(1,2);
    areaSolucao_candidata(9)=Solucao_candidata(1,2);
    areaSolucao_candidata(10)=Solucao_candidata(1,2);
    areaSolucao_candidata(11)=Solucao_candidata(1,2);
    areaSolucao_candidata(12)=Solucao_candidata(1,2);
    areaSolucao_candidata(13)=Solucao_candidata(1,3);
    areaSolucao_candidata(14)=Solucao_candidata(1,3);
    areaSolucao_candidata(15)=Solucao_candidata(1,3);
    areaSolucao_candidata(16)=Solucao_candidata(1,3);
    areaSolucao_candidata(17)=Solucao_candidata(1,3);
    areaSolucao_candidata(18)=Solucao_candidata(1,3);
    areaSolucao_candidata(19)=Solucao_candidata(1,3);
    areaSolucao_candidata(20)=Solucao_candidata(1,3);
    areaSolucao_candidata(21)=Solucao_candidata(1,3);
    areaSolucao_candidata(22)=Solucao_candidata(1,3);
    areaSolucao_candidata(23)=Solucao_candidata(1,3);
    areaSolucao_candidata(24)=Solucao_candidata(1,3); 

    
    ArqInput_C; %Arquivo com dados da Treliça
        
    index = find(areaSolucao_candidata>=areaCorte);
    areaAux=areaSolucao_candidata(index);
    conectividades=conectividades(index,:);    
    nBarrasAux=size(areaAux,2);
    
    areas = strcat('areas');
    save(areas,'Solucao_candidata');
    
    % Calculando deslocamentos não-lineares
    DeslVE = main();
    for id = 1:length(DeslVE)
        if mod(id,3) == 0
            z = DeslVE(id);
            DeslVE(id) = DeslVE(id-1);
            DeslVE(id-1) = z;
        end
    end
    
    % Adicionando às coordenadas nodais
    CoordVetor = CoordVetor + DeslVE';
    j = 1;
    for jd = 1:length(CoordVetor)
        if mod(jd,3) == 0
            CoordAux(j,1) = CoordVetor(jd-2);
            CoordAux(j,2) = CoordVetor(jd-1);
            CoordAux(j,3) = CoordVetor(jd);
            j = j+1;
        end
    end

    [nBarrasAux,areaAux,conectVetorAux,elemJuntaram,nodesC,numConst,conectividadesAux,elasticity,aalfa] = ROM(nBarrasAux,numNodes,CoordAux,conectividades,areaAux',nodesC); 
    nBarrasAux=size(areaAux,1);     

    [Desl,TensVE,Ke, Kg,Mmat]=Trelica_3d(nBarrasAux,areaAux,numNodes,numColConect,conectVetorAux,CoordVetor,nodesC,nodesL,numConst,elasticity,nodesCC,aalfa,ro,tamNodesCC);  %Calcula Deslocamentos e Tensões  

%      if max(abs(DeslVE)==99) || max(abs(DeslVE)==9999)     
%          lambdacrit=0;
%          ro=9999999;
%          freqnat=zeros(numNodes*2);
%      end


    lambda=sort(eig(Ke,-Kg));
    lambdapos=lambda(find(lambda>0));
    
    if isempty(lambdapos)==1
       lambdapos=0;
    end
    
    % Condição de distância
    ind1 = 0;
    ind2 = 1;
    lambdacrit(1) = lambdapos(1);
    
    if size(lambdapos,1)>1
    for ind1 = 2:size(lambdapos,1)
        if (abs(lambdapos(ind1) - lambdapos(ind1-1))) > 0.01
            ind2 = ind2 + 1; 
            lambdacrit(ind2,1) = lambdapos(ind1);
        end
    end
    end    

    if isempty(lambdacrit)==1
        lambdacrit=0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%

    frequencia=eig(Ke,Mmat);
    freq=sort(sqrt(frequencia)/(2*pi));
    
    % Condição de distância frequências
    ind1 = 0;
    ind2 = 1;
    freqnat(1) = freq(1);
    
    for ind1 = 2:size(freq,1)
        if (abs(freq(ind1) - freq(ind1-1))) > 0.01
            ind2 = ind2 + 1; 
            freqnat(ind2,1) = freq(ind1);
        end
    end
    
%     %%%%%%%%%%%%%%%%%%%%%
 
 
    clear Ke Kg Mmat; %TensAux;

    areaSolucao_candidata=areaAux';
     for j=1:nBarrasAux     

            x1=CoordAux(conectividadesAux(j,1),1);
            x2=CoordAux(conectividadesAux(j,2),1);
            y1=CoordAux(conectividadesAux(j,1),3);
            y2=CoordAux(conectividadesAux(j,2),3);
            z1=CoordAux(conectividadesAux(j,1),2);
            z2=CoordAux(conectividadesAux(j,2),2);

            L(j)=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);  

     end
    
    
    forca = 50000;
    f=funcpop(areaSolucao_candidata',ro,L,NumFO,lambdacrit(1),freqnat(1),forca,abs(DeslVE(2)));
%     f=funcpop(areaSolucao_candidata',ro,L,NumFO,lambdacrit(1:3),freqnat(1:3),forca,abs(DeslVE(2)));
    g=restricao(abs(DeslVE(2)),QuantRestr);


        
end
         

        
