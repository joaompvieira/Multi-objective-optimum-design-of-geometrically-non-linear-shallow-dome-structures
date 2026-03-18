function [f,g] = t60(Solucao_candidata)    

    Solucao_candidata=Solucao_candidata';
    dados_do_problema;
    ArqInput_C; %Arquivo com dados da Treliça
    
%     TensAux=zeros(nAreas,1);
    
%     tenZero=[1];
%     while isempty(tenZero)~=1
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
    areaSolucao_candidata(25)=Solucao_candidata(1,4);
    areaSolucao_candidata(26)=Solucao_candidata(1,4);
    areaSolucao_candidata(27)=Solucao_candidata(1,4);
    areaSolucao_candidata(28)=Solucao_candidata(1,4);
    areaSolucao_candidata(29)=Solucao_candidata(1,4);
    areaSolucao_candidata(30)=Solucao_candidata(1,4);
    areaSolucao_candidata(31)=Solucao_candidata(1,5);        
    areaSolucao_candidata(32)=Solucao_candidata(1,5);
    areaSolucao_candidata(33)=Solucao_candidata(1,5);
    areaSolucao_candidata(34)=Solucao_candidata(1,5);
    areaSolucao_candidata(35)=Solucao_candidata(1,5);
    areaSolucao_candidata(36)=Solucao_candidata(1,5);
    areaSolucao_candidata(37)=Solucao_candidata(1,5);
    areaSolucao_candidata(38)=Solucao_candidata(1,5);
    areaSolucao_candidata(39)=Solucao_candidata(1,5);
    areaSolucao_candidata(40)=Solucao_candidata(1,5);
    areaSolucao_candidata(41)=Solucao_candidata(1,5);
    areaSolucao_candidata(42)=Solucao_candidata(1,5);
    areaSolucao_candidata(43)=Solucao_candidata(1,6);
    areaSolucao_candidata(44)=Solucao_candidata(1,6);
    areaSolucao_candidata(45)=Solucao_candidata(1,6);
    areaSolucao_candidata(46)=Solucao_candidata(1,6);
    areaSolucao_candidata(47)=Solucao_candidata(1,6);
    areaSolucao_candidata(48)=Solucao_candidata(1,6);
    areaSolucao_candidata(49)=Solucao_candidata(1,7);
    areaSolucao_candidata(50)=Solucao_candidata(1,7);
    areaSolucao_candidata(51)=Solucao_candidata(1,7);
    areaSolucao_candidata(52)=Solucao_candidata(1,7);
    areaSolucao_candidata(53)=Solucao_candidata(1,7);
    areaSolucao_candidata(54)=Solucao_candidata(1,7);
    areaSolucao_candidata(55)=Solucao_candidata(1,7);
    areaSolucao_candidata(56)=Solucao_candidata(1,7);
    areaSolucao_candidata(57)=Solucao_candidata(1,7);
    areaSolucao_candidata(58)=Solucao_candidata(1,7);
    areaSolucao_candidata(59)=Solucao_candidata(1,7);
    areaSolucao_candidata(60)=Solucao_candidata(1,7);         
  

        index = find(areaSolucao_candidata>=areaCorte);
        areaAux=areaSolucao_candidata(index);
        conectividades=conectividades(index,:);    
        nBarrasAux=size(areaAux,2);
        
    areas = strcat('areas');
    save(areas,'Solucao_candidata');
        
    % Calculando deslocamentos não-lineares
    DeslV = main();
    for id = 1:length(DeslV)
        if mod(id,3) == 0
            z = DeslV(id);
            DeslV(id) = DeslV(id-1);
            DeslV(id-1) = z;
        end
    end
    
    % Adaptação dos nós
    DeslVE = DeslV;
        DeslVE(40:42) = 0;
        DeslVE(43:45) = DeslV(40:42);
        DeslVE(46:48) = 0;
        DeslVE(49:51) = DeslV(43:45);
        DeslVE(52:54) = 0;
        DeslVE(55:57) = DeslV(46:48);
        DeslVE(58:60) = 0;    
        DeslVE(61:63) = DeslV(49:51);
        DeslVE(64:66) = 0;
        DeslVE(67:69) = DeslV(52:54);
        DeslVE(70:72) = 0;
        DeslVE(73:75) = DeslV(55:57);

    %%%%%%%%%%%%
    
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
    
    %%%%%%%%%%%%%%%%%%%%%% 
    

    if isempty(lambdacrit)==1
        lambdacrit=0;
    end

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
    
    %%%%%%%%%%%%%%%%%%%%%%

        
%         tenZero = find(TensVE==0)'; % all zeros
%  
% 
%         if isempty(tenZero)~=1   
%             for j=1:size(tenZero,1)
%                 [l,r]=find(conectividades(:,1)==conectividadesAux(tenZero(j),1) & conectividades(:,2)==conectividadesAux(tenZero(j),2));
%                 if isempty(l)==1
%                     [l,r]=find(   all(round((TensVE(:,size(TensVE,2)-size(elemJuntaram,2)+1:size(TensVE,2))))==0)    );
%                     for k=1:size(r,2)
%                         aux=index(elemJuntaram{r(k)});
%                         for kk=1:size(aux,2)  
%                             Solucao_candidata(1,aux(kk))=areaMin;
%                         end
%                     end
%                 else                    
%                     for k=1:size(r,2)
%                         aux=index(l(k));  
%                         Solucao_candidata(1,aux)=areaMin;
%                     end         
%                 end
%             end      
%         end
%         
%         if isempty(elemJuntaram)~=1
%             for jj=1:size(elemJuntaram,2)
%                 TensAuxx=TensVE(size(TensVE,1)-size(elemJuntaram{jj},1)+1:size(TensVE,1));
%                 TensVE=TensVE(1:size(TensVE,1)-size(elemJuntaram{1},1));
% 
%                 for j=1:size(elemJuntaram{jj},1)        
%                     aux=index(setdiff(elemJuntaram{jj}(j,:),0,'stable'));
%                     index=setdiff(index,index([setdiff(elemJuntaram{jj}(j,:),0,'stable')]),'stable');
%                     index=[index aux];                      
%                     if isempty(TensVE)~=1
%                         TensVE(size(TensVE,1)+1:size(TensVE,1)+size(aux,2))=TensAuxx;
%                     else
%                         TensVE(size(TensVE,1):size(aux,2))=TensAuxx;
%                     end
%                     clear aux;
%                 end
%             end
%             TensAux(index)=TensVE;
%             TensVE=TensAux;
%         else
%             TensAux(index)=TensVE;
%             TensVE=TensAux;
%         end
%             
%         
%         
%     end            


    clear Ke Kg Mmat;% TensAux;

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
    
    
    
    forca = 10000;
    f=funcpop(areaSolucao_candidata',ro,L,NumFO,lambdacrit(1),freqnat(1),forca,abs(DeslVE(2)));
%     f=funcpop(areaSolucao_candidata',ro,L,NumFO,lambdacrit(1),freqnat(1));
%     f=funcpop(areaSolucao_candidata',ro,L,NumFO);
    g=restricao(abs(DeslVE),QuantRestr);


        
end
         

        
