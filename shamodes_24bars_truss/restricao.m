function v=restricao(desl,QuantRestr);    
    v=zeros(1,QuantRestr); 

    
%     for j=1:1:nBarras  
%         v(j) =  (abs(Tens(j))/2.758e+8) - 1;     %40ksi   
%     end    
    
v(1)=(desl/0.67)-1; 

%     v(32)=1.0-(lambdacrit/1.0);
%     v(32)=1.0-(freqnat/10.0);

    % Se for menor que zero é porque não viola as restrições e portanto recebe
    % violação = 0

    for j=1:QuantRestr
        if v(j)<=0          
            v(j)=0;
        end    
    end    

end
