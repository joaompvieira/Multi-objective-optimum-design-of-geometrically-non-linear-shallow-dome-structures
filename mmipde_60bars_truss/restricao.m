function v=restricao(desl,QuantRestr);    
    v=zeros(1,QuantRestr); 


    v(1) = (desl(2)/0.25)-1; % nó 1 dir z
    v(2) = (desl(5)/0.25)-1; % nó 2 dir z
    v(3) = (desl(23)/0.25)-1; % nó 8 dir z
    v(4) = (desl(44)/0.25)-1; % nó 15 dir z


    % Se for menor que zero é porque não viola as restrições e portanto recebe
    % violação = 0

    for j=1:QuantRestr
        if v(j)<=0          
            v(j)=0;
        end    
    end    

end
