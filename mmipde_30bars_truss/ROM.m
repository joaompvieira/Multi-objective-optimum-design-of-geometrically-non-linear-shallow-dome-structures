function [nBarrasAux,areaAux,conectVetorAux,elemJuntaram,nodesC,numConst,conectividadesAux,elasticity,aalfa] = ROM (nBarras,numNodes,CoordAux,conectividades,area,nodesC)

% area=[0.000100000000000000;0;0.000300000000000000;0;0.000150000000000000;0.000150000000000000;0;0;0;0.000150000000000000;0;0;0;0;0;0.000100000000000000;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.000100000000000000;0;0;0;0;0;0;0;0;0;0;0;0.000100000000000000;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.00200000000000000;0.000150000000000000;0.00185000000000000;0.00190000000000000;0;0.000200000000000000;0;0;0.000150000000000000;0;0;0;0;0;0.000150000000000000;0.000100000000000000;0.000350000000000000;0.000250000000000000;0;0.000100000000000000;0;0;0;0;0;0;0;0.000100000000000000;0;0.000150000000000000;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.000550000000000000;0;0;0;0;0;0;0;0;0;0.00190000000000000;0.000700000000000000;0;0;0;0.000100000000000000;0;0;0;0;0;0.000300000000000000;0;0.00195000000000000;0;0;0.000100000000000000;0;0;0;0;0.000150000000000000;0;0;0.000450000000000000;0;0.00145000000000000;0;0;0;0.000150000000000000;0;0;0;0;0;0;0;0.000700000000000000;0;0;0;0;0;0;0.000250000000000000;0.000700000000000000;0;0.000450000000000000;0;0.000100000000000000;0;0;0;0;0;0.00195000000000000;0;0;0.000200000000000000;0;0;0;0;0;0;0;0.00145000000000000;0;0;0;0;0;0;0;0;0;0.000700000000000000;0;0;0;0;0.00190000000000000;0.000150000000000000;0.000100000000000000;0;0;0;0;0.000300000000000000;0.00185000000000000;0.000350000000000000;0;0;0.000150000000000000;0;0;0.00190000000000000;0.000250000000000000;0;0;0;0;0;0;0;0.000550000000000000;0;0;0.00200000000000000;0;0;0.000100000000000000;0;0.000150000000000000;0.000300000000000000;0;0;0;0.000150000000000000;0;0;0;0;0.000150000000000000;0;0;0;0.000100000000000000;0.000100000000000000;0];

[l,r] = find(area~=0);


nosQueSaem2=[];
conectAux=0;
conectAux=[conectividades(l,1) conectividades(l,2)];
areaAux=0;
areaAux=area(find(area~=0));


% for j=1:1:size(areaAux,1)
%     linha=line([CoordAux(conectividades(l(j),1),1) CoordAux(conectividades(l(j),2),1)],[CoordAux(conectividades(l(j),1),3) CoordAux(conectividades(l(j),2),3)],[CoordAux(conectividades(l(j),1),2) CoordAux(conectividades(l(j),2),2)]);            
% end
% hold on;
% 
% raioesfera=6;
% b=unique(conectividades(l,1:2));
% for j=1:size(b,1)
%     text(CoordAux(b(j),1)+0.2,CoordAux(b(j),3)+0.2,CoordAux(b(j),2)+0.1,num2str(b(j)));
%     [X2,Y2,Z2]=sphere;
%      X2=X2/raioesfera;
%      Y2=Y2/raioesfera;
%      Z2=Z2/raioesfera;
%      a=surf(X2+CoordAux(b(j),1),Y2+CoordAux(b(j),3),Z2+CoordAux(b(j),2));
%      a.MeshStyle='row';
%     a.LineStyle='non';
%     a.FaceColor=[0 0 0];
%     a.FaceAlpha=1;
%     hold on;
% end
% % keyboard;

cont=1;
cont2=1;
for i=1:numNodes
    [l,r] = find(conectAux(:,1)==i | conectAux(:,2)==i);
%     if isempty(l)==1
%         nosQueSaem(cont)=i;
%         cont=cont+1;    
    if isempty(l)==0
        tan=zeros(1,size(l,1));
        for j=1:size(l,1)
            [elem,col]=find(((conectividades(:,1)==conectAux(l(j),1))&(conectividades(:,2)==conectAux(l(j),2))) | ((conectividades(:,2)==conectAux(l(j),1))&(conectividades(:,1)==conectAux(l(j),2))));
            x1=CoordAux(conectividades(elem,1),1);
            x2=CoordAux(conectividades(elem,2),1);
            y1=CoordAux(conectividades(elem,1),3);
            y2=CoordAux(conectividades(elem,2),3);
            tan(j)=(y2-y1)/(x2-x1);            
        end
        tan=unique(tan);
        if size(tan,2)==1 & size(l,1)>1
%             area1=areaAux(l(1,1));
%             area2=areaAux(l(2,1));
%             if (area1 == area2)
%                 keyboard;
                nosQueSaem2(cont2)=i;
                cont2=cont2+1;
%             end
        end        
    end 
        
end

elemJuntaram=[];
% nosQueSaem
% nosQueSaem2
nosTan2=[];
if isempty(nosQueSaem2)==0
    clear tan;
    [l,r] = find(area~=0);
    for j=1:nBarras
        if area(j)~=0
            x1=CoordAux(conectividades(j,1),1);
            x2=CoordAux(conectividades(j,2),1);
            y1=CoordAux(conectividades(j,1),3);
            y2=CoordAux(conectividades(j,2),3);
            tan(j)=(y2-y1)/(x2-x1);        
        else
            tan(j)=NaN;
        end
    end


    nosQueFicam=setdiff([1:numNodes],nosQueSaem2);
    % nosQueSaem2=setdiff(nosQueSaem2,[1 2]);

    tanAux=unique(tan);
    nosTan=[];
    l=isnan(tanAux);
    l=find(l==0);

    tanAux=tanAux(l)';

    cont=1;
    for i=1:size(tanAux,1)
        elemTan=find(tan==tanAux(i));    
        if size(elemTan,2)>1
            elemTanDefinitivo{cont}=[];
            for j=1:size(nosQueSaem2,2)
                elemTanDefinitivo{cont}=[unique(elemTanDefinitivo{cont}) elemTan(find(conectividades(elemTan,1)==nosQueSaem2(j) | conectividades(elemTan,2)==nosQueSaem2(j)))];
            end
            
            if isempty(elemTanDefinitivo{cont})==0            
                nosTan{cont}=unique(conectividades(elemTanDefinitivo{cont},:));
                cont=cont+1;
                
            end
        end
    end
    
    cont=1;
    for kkk=1:size(nosTan,2)            %%%%%%%%%%% SEPARANDO CORRENTES COM MESMA TANGENTE  %%%%%%%%%%%%%%%%%%%%%%%% 
        for j=2:size(nosTan{kkk},1)
            x1=CoordAux(nosTan{kkk}(1),1);
            x2=CoordAux(nosTan{kkk}(j),1);
            y1=CoordAux(nosTan{kkk}(1),3);
            y2=CoordAux(nosTan{kkk}(j),3);
            tanNo(j-1)=((y2-y1)/(x2-x1));                    
        end              
%         auxTan=conectAux([find(conectAux(:,1)==nosTan{kkk}(1) | conectAux(:,2)==nosTan{kkk}(1))],:);
%         for jj=1:size(auxTan,1)
%             for jjj=1:size(nosQueSaem2,2)
%                 if (conectAux(auxTan(jj),1)==nosQueSaem2(jjj) | conectAux(auxTan(jj),2)==nosQueSaem2(jjj));
%                     auxTan=conectAux(auxTan(jj),:);
%                     break;
%                 end
%             end
%              
%         end
        
        tanAlinhamento=unique(tanNo); 

        while size(tanAlinhamento,2)>1
            for j=1:size(tanAlinhamento,2)
                noTanQqer=[];
                for jj=1:size(nosQueSaem2,2)
                    for jjj=1:size(nosTan{kkk},1)
                        if nosTan{kkk}(jjj)==nosQueSaem2(jj)
                            noTanQqer=nosQueSaem2(jj);
                        end
                    end
                end
                [l2,r2]=find(conectAux(:,1)==noTanQqer | conectAux(:,2)==noTanQqer);
                x1=CoordAux(conectAux(l2(1),1),1);
                x2=CoordAux(conectAux(l2(1),2),1);
                y1=CoordAux(conectAux(l2(1),1),3);
                y2=CoordAux(conectAux(l2(1),2),3);
                tanQqerAlinhamento=((y2-y1)/(x2-x1));
                
                [l,r]=find(tanNo==tanQqerAlinhamento);
                flag=0;
%                 for jj=1:size(conectAux,1)
%                     if ((conectAux(jj,1)==nosTan{kkk}(1) & conectAux(jj,2)==nosTan{kkk}(r(1)+1)) | (conectAux(jj,2)==nosTan{kkk}(1) & conectAux(jj,1)==nosTan{kkk}(r(1)+1)))
%                         flag=1;
%                         break;
%                     end
%                 end
                if size(r,2)>1 
                    nosTan{size(nosTan,2)+1}=nosTan{kkk}([1 r+1]);
                    nosTan{kkk}=nosTan{kkk}(setdiff(1:size(nosTan{kkk},1),[1 r+1]));  
                    break;
                end
                clear noTanQqer tanQqerAlinhamento l2 r2;
            end                    
            clear tanNo tanAlinhamento;
            for j=2:size(nosTan{kkk},1)
                x1=CoordAux(nosTan{kkk}(1),1);
                x2=CoordAux(nosTan{kkk}(j),1);
                y1=CoordAux(nosTan{kkk}(1),3);
                y2=CoordAux(nosTan{kkk}(j),3);
                tanNo(j-1)=((y2-y1)/(x2-x1));                    
            end 
            tanAlinhamento=unique(tanNo);                   

        end
        clear tanNo tanAlinhamento;
    end

        %%%%%%%% FIM DE SEPARANDO AS CORRENTES COM MESMA TANGENTE  %%%%%%%%%% 

        %%%%%%% JUNTAR MAIS DE UM NO NO MESMO ALINHAMENTO %%%%%%%%%%%%%        

    for j=1:size(nosTan,2)  

        nosTanAux=[];
        cont3=1;
%             nosTan{j}=[];
        cont2=1;

        for jj=1:size(nosTan{j},1)
            for jjj=1:size(nosQueSaem2,2)
                if  nosTan{j}(jj)==nosQueSaem2(jjj)
                    nosSaemAlinhamento(cont2)=nosTan{j}(jj);
                    cont2=cont2+1;
                    break;
                end
            end
        end                       
        while isempty(nosSaemAlinhamento)~=1
            [l,r]=find(conectAux(:,1)==nosSaemAlinhamento(1) | conectAux(:,2)==nosSaemAlinhamento(1));
            nosSaemAlinhamento=setdiff(nosSaemAlinhamento,nosSaemAlinhamento(1));
            aux=unique(conectAux(l,:));
            cont2=0;
            for jj=1:size(aux,1)
                if (isempty(find(nosSaemAlinhamento==aux(jj)))~=1)
                    cont2=cont2+1;
                end
            end
            while cont2>=1
                [l,r]=find(conectAux(:,1)==nosSaemAlinhamento(1) | conectAux(:,2)==nosSaemAlinhamento(1));
                nosSaemAlinhamento=setdiff(nosSaemAlinhamento,nosSaemAlinhamento(1));
                aux=[aux;unique(conectAux(l,:))];
                cont2=0;
                for jj=1:size(aux,1)
                    if (isempty(find(nosSaemAlinhamento==aux(jj)))~=1)
                        cont2=cont2+1;
                    end
                end
            end           
        nosTanAux{cont3}=unique(aux);
        cont3=cont3+1;
        end    

        nosTan2=[nosTan2 nosTanAux]; 
    end  
    cont=cont+1;
end

clear tan;
clear tanAux;

if (isempty(nosTan2)==0)
    nosTan=nosTan2;
    clear nosTan2;



    elemJuntaram=[];

%         clear nosTan;
    cont=1;
    for i=1:size(nosTan,2)
        aux=unique(nosTan{i})';           
        
        elementos=[];
        for j=1:size(aux,2)
            for jj=j+1:size(aux,2)
                elementos=[elementos;find(conectAux(:,1)==aux(j) & conectAux(:,2)==aux(jj))];
            end
        end
        elementos=conectAux(elementos,:);
        nos=unique(elementos);
        cont2=zeros(size(nos,1),1);
        for j=1:size(nos,1)
            for jj=1:size(elementos,1)
                for jjj=1:size(elementos,2)
                    if (elementos(jj,jjj)==nos(j))
                        cont2(j)=cont2(j)+1; 
                    end
                end
            end
        end
        [l,r]=find(cont2==1);
        noMin=l(1);
        noMax=l(2);
        nosTan{i}=[aux(noMin(1,1)) aux(noMax(1,1))];
        conectAux=[conectAux;[nosTan{i}]];
        aux=setdiff(aux,[aux(noMin) aux(noMax)]);

        for j=1:size(aux,2)        

            l = find(conectAux(:,1)==aux(j) | conectAux(:,2)==aux(j));
            for k=1:size(l,1)
                l2(k)=find(conectividades(:,1)==conectAux(l(k),1) & conectividades(:,2)==conectAux(l(k),2));
                elemJuntaram(i,cont)=l2(k);
                cont=cont+1;
            end
            conectAux=conectAux(setdiff([1:size(conectAux,1)]',[l]),:);
            areaAux=areaAux(setdiff([1:size(areaAux,1)]',[l]),:);
            a1=min(area(l2));
        end


        areaAux=[areaAux;a1];
        

    end
end



if isempty(elemJuntaram)~=1
    for j=1:size(elemJuntaram,1)
        elemJuntaram2{j}=setdiff(elemJuntaram(j,:),0);
    end
    elemJuntaram=elemJuntaram2;
    clear elemJuntaram2;
end

conectividades=conectAux;
nBarrasAux=size(conectAux,1);
conectividadesAux=conectividades;
  
conectVetorAux=[];
for cLncv=1:size(conectividades,1)
    conectVetorAux=[conectVetorAux conectividades(cLncv,:)];
end  

nosSaem=setdiff([1:numNodes],unique(unique(conectividades)));
for j=1:size(nosSaem,2)
    nodesC(3*nosSaem(j)-2:3*nosSaem(j))=[1 1 1];
end

numConst=size(find(nodesC==1),2);


elasticity=[];
aalfa=[];
for aaaaux=1:nBarrasAux
    elasticity(1,aaaaux)=10000000;
    
end                 

aalfa=zeros(1,nBarrasAux);
