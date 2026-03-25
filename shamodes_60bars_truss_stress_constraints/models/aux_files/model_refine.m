%% model_refine.m
% Auxiliary script to refine elements discretization of a given model.
%
% Select the input file to be refined (variable 'file') and the number of
% subdivisions of each element (variable 'subdiv').
%
% The input file of the original model will not be not edited.
% This script generates a new input file with the data of the refined model.
% The new file has the same name of the original file plus the sufix
% '_refined'.
%
% The number IDs of the nodes of the original model remain the same.
% The number IDs of the new nodes, added to the elements to refine the model,
% start from the last number ID of the original model.
%
% Each subdivision of the refined model is a new element and, therefore,
% the number IDs of the elements are different from the original model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the new refined input file, element loads (distributed forces or 
% temperature gradient) are reseted, i.e., if the original model has element
% loads, they must be added again to the refined model, considering the new
% element ID numbers!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = 'path_name/file_name';
subdiv = 10;

%% Main code
clc; close all; clearvars -except file subdiv

% Check for a valid number of subdivision
if (subdiv <= 1 || ~(floor(subdiv) == subdiv))
    fprintf('invalid number of subdivisions!\nIt must be an integer greater than 1...\n');
    return;
end

% Open input file
fin = fopen(file,'rt');
if (fin < 0)
    fprintf('Error opening input file!\n');
    return;
end

% Get node coordinates and element connectivity
[status,coords] = getCoords(fin);
if (status == 0)
    fclose(fin);
    return;
end
[status,connect] = getConect(fin,size(coords,1));
if (status == 0)
    fclose(fin);
    return;
end

% Refine elements
[newCoords,newConnect] = refineModel(subdiv,coords,connect);

% Open output file
namein  = fopen(fin);
nameout = strcat(namein(1:end-4),'_refined.txt');
fout    = fopen(nameout,'w');

% Write output file
createRefinedFile(fin,fout,subdiv,newCoords,newConnect);

% Draw new model
drawRefinedModel(newCoords,newConnect);

% Close files
fclose(fin);
fclose(fout);

%% Auxiliary Functions
%--------------------------------------------------------------------------
function [status,coords] = getCoords(fin)
status = 0;
while ~feof(fin)
    string = deblank(fgetl(fin));
    if (strcmp(string,'%MODEL.NODE_COORD'))
        [status,coords] = readNodeCoord(fin);
        return;
    end
end
end

%--------------------------------------------------------------------------
function [status,coords] = readNodeCoord(fin)
status = 1;
nnp = fscanf(fin,'%d',1);
if (nnp <= 0)
    fprintf('Invalid number of node coordinates!\n');
    status = 0;
    return;
end
coords = zeros(nnp,3);
for i = 1:nnp
    id = fscanf(fin,'%d',1);
    if (id <= 0 || id > nnp)
        fprintf('Invalid node ID for coordinate specification!\n');
        status = 0;
        return;
    end
    [values,count] = fscanf(fin,'%f',3);
    if (count ~= 3)
        fprintf('Invalid node coordinate!\n');
        status = 0;
        return;
    end
    coords(id,:) = values;
end
end

%--------------------------------------------------------------------------
function [status,connect] = getConect(fin,nnp)
status = 0;
while ~feof(fin)
    string = deblank(fgetl(fin));
    if (strcmp(string,'%MODEL.ELEMENT'))
        [status,connect] = readConnect(fin,nnp);
        return;
    end
end
end

%--------------------------------------------------------------------------
function [status,connect] = readConnect(fin,nnp)
status = 1;
nel = fscanf(fin,'%d',1);
if (nel <= 0)
    fprintf('Invalid number of elements!\n');
    status = 0;
    connect = [];
    return;
end
connect = zeros(nel,2);
for i = 1:nel
    id = fscanf(fin,'%d',1);
    if (id <= 0 || id > nel)
        fprintf('Invalid element ID!\n');
        status = 0;
        return;
    end
    [values,count] = fscanf(fin,'%f',10);
    if (count ~= 10 ||...
        (values(4) <= 0 && values(4) > nnp)  ||...
        (values(5) <= 0 && values(5) > nnp))
        fprintf('Invalid element properties!\n');
        status = 0;
        return;
    end
    connect(id,:) = [values(4),values(5)];
end
end

%--------------------------------------------------------------------------
function [newCoords,newConnect] = refineModel(subdiv,oldCoords,oldConnect)
old_nnp = size(oldCoords,1);
old_nel = size(oldConnect,1);
new_nnp = old_nnp + old_nel * (subdiv - 1);
new_nel = old_nel * subdiv;
newCoords  = zeros(new_nnp,3);
newConnect = zeros(new_nel,2);
newCoords(1:old_nnp,:) = oldCoords;
n = old_nnp + 1; % new nodes iterator
e = 1;           % new elements iterator
for i = 1:old_nel
    x1 = oldCoords(oldConnect(i,1),1);
    y1 = oldCoords(oldConnect(i,1),2);
    z1 = oldCoords(oldConnect(i,1),3);
    x2 = oldCoords(oldConnect(i,2),1);
    y2 = oldCoords(oldConnect(i,2),2);
    z2 = oldCoords(oldConnect(i,2),3);
    dx = (x2 - x1) / subdiv;
    dy = (y2 - y1) / subdiv;
    dz = (z2 - z1) / subdiv;
    for j = 1:subdiv-1
        x = x1 + j*dx;
        y = y1 + j*dy;
        z = z1 + j*dz;
        newCoords(n,:) = [x,y,z];
        if (j == 1)
            newConnect(e,:) = [oldConnect(i,1),n];
        else
            newConnect(e,:) = [n-1,n];
        end
        n = n + 1;
        e = e + 1;
    end
    newConnect(e,:) = [n-1,oldConnect(i,2)];
    e = e + 1;
end
end

%--------------------------------------------------------------------------
function createRefinedFile(fin,fout,subdiv,coords,connect)
frewind(fin)
while ~feof(fin)
    line = fgets(fin);
    fwrite(fout, line);
    if (strcmp(deblank(line),'%MODEL.NODE_COORD'))
        oldNodes = fscanf(fin,'%d',1);
        break;
    end
end
fprintf(fout,'%d\n',size(coords,1));
for i = 1:size(coords,1)
    fprintf(fout,'%d   %f %f %f\n',i,coords(i,1),coords(i,2),coords(i,3));
end
for i = 1:oldNodes+1
    fgets(fin);
end
while ~feof(fin)
    line = fgets(fin);
    fwrite(fout, line);
    if (strcmp(deblank(line),'%MODEL.ELEMENT'))
        oldElems = fscanf(fin,'%d',1);
        break;
    end
end
m = 1; % new elements iterator
fprintf(fout,'%d\n',size(connect,1));
for i = 1:oldElems
    props = fscanf(fin,'%d',11);
    for j = 1:subdiv
        % ID
        fprintf(fout,'%d     ',m);
        % Type
        fprintf(fout,'%d   ',props(2));
        % Material /Section
        fprintf(fout,'%d ',  props(3));
        fprintf(fout,'%d   ',props(4));
        % Nodes
        fprintf(fout,'%d ',  connect(m,1));
        fprintf(fout,'%d   ',connect(m,2));
        % End liberations
        if (j == 1)
            fprintf(fout,'%d ',  props(7));
            fprintf(fout,'%d   ',1);
        elseif (j == subdiv)
            fprintf(fout,'%d ',  1);
            fprintf(fout,'%d   ',props(8));
        else
            fprintf(fout,'%d ',  1);
            fprintf(fout,'%d   ',1);
        end
        % Orientation vector
        fprintf(fout,'%d ',props(9));
        fprintf(fout,'%d ',props(10));
        fprintf(fout,'%d\n', props(11));
        m = m + 1;
    end
end
oldLoads = 0;
while ~feof(fin)
    line = fgets(fin);
    fwrite(fout, line);
    if (strcmp(deblank(line),'MODEL.LOAD.ELEM_UNIFORM'))
        oldLoads = fscanf(fin,'%d',1);
        break;
    end
end
for i = 1:oldLoads+1
    fgets(fin);
end
oldLoads = 0;
while ~feof(fin)
    line = fgets(fin);
    fwrite(fout, line);
    if (strcmp(deblank(line),'MODEL.LOAD.ELEM_LINEAR'))
        oldLoads = fscanf(fin,'%d',1);
        break;
    end
end
for i = 1:oldLoads+1
    fgets(fin);
end
while ~feof(fin)
    line = fgets(fin);
    fwrite(fout, line);
end
end

%--------------------------------------------------------------------------
function drawRefinedModel(coords,connect)
figure('Name','refined_model');
hold on
grid off
axis equal
title('Refined Model');
for i = 1:size(connect,1)
    x1 = coords(connect(i,1),1);
    y1 = coords(connect(i,1),2);
    z1 = coords(connect(i,1),3);
    x2 = coords(connect(i,2),1);
    y2 = coords(connect(i,2),2);
    z2 = coords(connect(i,2),3);
    x = [x1,x2];
    y = [y1,y2];
    z = [z1,z2];
    scatter3(x,y,z,10,'filled','black');
    plot3(x,y,z,'black');
end
xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
zLimits = get(gca,'ZLim');
gapX    = (xLimits(2)-xLimits(1))/10;
gapY    = (yLimits(2)-yLimits(1))/10;
gapZ    = (zLimits(2)-zLimits(1))/10;
minX    = xLimits(1) - gapX;
maxX    = xLimits(2) + gapX;
minY    = yLimits(1) - gapY;
maxY    = yLimits(2) + gapY;
minZ    = zLimits(1) - gapZ;
maxZ    = zLimits(2) + gapZ;
xlim([minX maxX]);
ylim([minY maxY]);
zlim([minZ maxZ]);
if (gapZ > 0)
    view(3);
end
end
