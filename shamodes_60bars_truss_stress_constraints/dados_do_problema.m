DimEspa=7; % Dimensão do Espaço
nAreas=7;
NumFO=4; % Número de Funções Objetivo

nBarras=60;

QuantRestr=64; % Quantidade de Restrições  %MOSO 3

LimSup(1:nAreas)=0.0020; % Limite Superior das Variáveis
LimInf(1:nAreas)=0.00005; % Limite Inferior das Variáveis

ProbCruzamento=0.9; % Probabilidade de Cruzamento da Evolução Diferencial
F=0.3; % Fator F da Evolução Diferencial
NumExec=10; % Número de Execuções Independentes
TamanhoPopulacao=20; % Tamanho da População
NumMaxGeracoes=100; % Numero máximo de gerações
areaCorte=0.00001;
areaMin=0;