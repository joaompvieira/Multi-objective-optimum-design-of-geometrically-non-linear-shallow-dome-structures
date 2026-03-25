function resultado=funcpop(A,ro,L,NumFO,lambdacrit,freqnat,carga,d)
    
 
%     resultado(1,1)=ro*(L*A);       % MOSOP1
    
 
%     resultado(1,1)=ro*(L*A);        % MOSOP2 
%     resultado(1,2)=-lambdacrit;
%     resultado(1,3)=-freqnat;

      resultado(1,1)=ro*(L*A);        % MOSOP3
      resultado(1,2)=-lambdacrit;
      resultado(1,3)=-freqnat;
      resultado(1,4)=carga*d;
   
   
     resultado=resultado(1,1:NumFO);

end
