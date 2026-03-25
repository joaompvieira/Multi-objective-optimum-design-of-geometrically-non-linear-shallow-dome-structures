% Copyright (c) 2019, Natee Panagant
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of  nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prob='MOSOP3_nlinear_stress';
nloop=100;
nsol=10;
narchive=2*nsol;
nvar=7;
nrun=10; 

lb = [0.00005;0.00005;0.00005;0.00005;0.00005;0.00005;0.00005]; 
ub = [0.0020;0.0020;0.0020;0.0020;0.0020;0.0020;0.0020];

% lb=zeros(nvar,1);
% ub=ones(nvar,1);
%algo={'SHAMODE'};
% algo={'SHAMODE_WO'};

 algo={'SHAMODE';
       'SHAMODE_WO'};
%prob = 'MOSO4';
fun='t60';
% algo={'SHAMODE'};
rst=cell(numel(algo),1);
%for i=1 % No. of Objective
% runningTime=zeros(2,1);
    for j=1:numel(algo) % No. of Algorithm
        algoj=algo{j};
        HV=[];
        tic;
        for k=1:nrun % No. of Run                    
            rst{j}=feval(algoj,fun,nloop,nsol,nvar,narchive,lb,ub,k);
            final_candidatas = rst{j}.ppareto{end}';
            final_objetivos = rst{j}.fpareto{end}';
            for i=1:size(final_objetivos,1)
                aux(i,:)=[final_candidatas(i,:) final_objetivos(i,:)];
            end
            if j==1
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 resultados=strcat('Resultados_',fun,'_',prob,'_SHAMODE.mat');
                Execucao_SHAMODE{k} = aux;
                resultados=('MOSOP3_nlinear_stress_SHAMODE');
                save(resultados,'Execucao_SHAMODE');
            else
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 resultados=strcat('Resultados_',fun,'_',prob,'_SHAMODE_WO.mat'); 
                 Execucao_SHAMODEWO{k} = aux;
                 resultados=('MOSOP3_nlinear_stress_SHAMODEWO');
                 save(resultados,'Execucao_SHAMODEWO');
            end
           % save(resultados,'Execucao');
        end
%         runningTime(j)= toc;

   end
%end