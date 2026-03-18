%% model_arch.m
% Auxiliary script to generate 2D coordinate points for a circular arch model.
% Select the span length (s), height (h), and number of elements to
% dicretize the arch (n), and the script will generate an output file
% entitled arch_model.txt with pairs of nodal coordinates.
% The script also draws the model.

s = 2.00;   % span length
h = 1.00;   % height 
n = 32;     % number of elements


%% Main code
clc; clearvars -except s h n

% Generate coordinates
x = zeros(n+1,1);
y = zeros(n+1,1);
r = (s^2/4+h^2)/(2*h);
a = atan2((s/2),(r-h));
ang = pi-(pi/2-a);
inc = 2*a/n;
for i=1:n+1
    x(i)=r*cos(ang);
    y(i)=r*sin(ang)-(r-h);
    ang=ang-inc;
end

% Plor arch
hold on
axis equal
for i = 1:n
    x1 = x(i);
    y1 = y(i);
    x2 = x(i+1);
    y2 = y(i+1);
    scatter([x1,x2],[y1,y2],10,'filled','black');
    plot(x,y,'black');
end

% Write coordinate vectors
fid = fopen('arch_model.txt','w');
if (fid < 0)
    fprintf('Error opening input file!\n');
    return;
end
for i=1:n+1
    fprintf(fid,'%d   %.15f   %.15f\n',i,x(i),y(i));
end
fclose(fid);
