function v=restricao(desl,Tens,nBarras,A,L,QuantRestr)    
    
    v=zeros(1,QuantRestr);
    k = 3.96; % Euler buckling coefficient
    E = 73e9; % Elasticity
    
    for j=1:nBarras 
        if Tens(j) >= 0
            v(j) =  (abs(Tens(j))/150e6) - 1; 
        else
            sig = (k*E*A(j))/((L(j))^2);
            if sig <= 150e6
               v(j) = (abs(Tens(j))/sig) - 1; 
            else
               v(j) =  (abs(Tens(j))/150e6) - 1;
            end
        end
    end    

    v(j+1) = (desl(2)/0.25)-1; % nó 1 dir z
    v(j+2) = (desl(5)/0.25)-1; % nó 2 dir z
    v(j+3) = (desl(23)/0.25)-1; % nó 8 dir z
    v(j+4) = (desl(44)/0.25)-1; % nó 15 dir z


    % Se for menor que zero é porque não viola as restrições e portanto recebe
    % violação = 0

    for i=1:QuantRestr
        if v(i)<=0          
            v(i)=0;
        end    
    end

end
