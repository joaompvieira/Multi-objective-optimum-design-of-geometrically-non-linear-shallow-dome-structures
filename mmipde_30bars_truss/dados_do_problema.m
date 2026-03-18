DimEspa=4; % Dimensão do Espaço
nAreas=4;
NumFO=4; % Número de Funções Objetivo

nBarras=30;

% QuantRestr=nBarras+1; % Quantidade de Restrições  %MOSO 1
% QuantRestr=nBarras+1; % Quantidade de Restrições  %MOSO 2
QuantRestr=2; % Quantidade de Restrições  %MOSO 3

LimSup(1:nAreas)=2; % Limite Superior das Variáveis
LimInf(1:nAreas)=0.1; % Limite Inferior das Variáveis
% LimSup(1,nAreas+1)=60*2.54/100; % Limite Superior das Variáveis
% LimInf(1,nAreas+1)=20*2.54/100; % Limite Inferior das Variáveis
% LimSup(1,nAreas+2)=80*2.54/100; % Limite Superior das Variáveis
% LimInf(1,nAreas+2)=40*2.54/100; % Limite Inferior das Variáveis
% LimSup(1,nAreas+3)=130*2.54/100; % Limite Superior das Variáveis
% LimInf(1,nAreas+3)=90*2.54/100; % Limite Inferior das Variáveis
% LimSup(1,nAreas+4)=80*2.54/100; % Limite Superior das Variáveis
% LimInf(1,nAreas+4)=40*2.54/100; % Limite Inferior das Variáveis
% LimSup(1,nAreas+5)=140*2.54/100; % Limite Superior das Variáveis
% LimInf(1,nAreas+5)=100*2.54/100; % Limite Inferior das Variáveis
% LimSup(1,nAreas+1:DimEspa)=1000*2.54/100*ones(1,3); % Limite Superior das Variáveis
% LimInf(1,nAreas+1:DimEspa)=180*2.54/100*ones(1,3); % Limite Inferior das Variáveis
ProbCruzamento=0.9; % Probabilidade de Cruzamento da Evolução Diferencial
F=0.3; % Fator F da Evolução Diferencial
NumExec=10; % Número de Execuções Independentes
TamanhoPopulacao=20; % Tamanho da População
NumMaxGeracoes=100; % Numero máximo de gerações
% filtro=0.2;
areaCorte=0.001;
areaMin=0;