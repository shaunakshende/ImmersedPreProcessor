function NodeSetGenerator(input_file)
fileID = fopen(input_file,'r');
formatSpec = '%f %f %f %d %f\n';
sizeCoor = [5 Inf];
Coor = fscanf(fileID, formatSpec, sizeCoor)';
fclose(fileID);
X = Coor(:,1); Y = Coor(:, 2); 
ID = 1:length(X);

Xmin = min(X);
Ymin = min(Y);
Xmax = max(X);
Ymax = max(Y);

nodeset_left = ID(X==Xmin)';
nodeset_right = ID(X==Xmax)';

nodeset_bottom = ID(Y==Ymin)';
nodeset_top = ID(Y==Ymax)';

% Save Left:
fileID = fopen('Left.txt','w');
fprintf(fileID,'%d\n',nodeset_left);
fclose(fileID);

% Save Right:
fileID = fopen('Right.txt','w');
fprintf(fileID,'%d\n',nodeset_right);
fclose(fileID);

% Save Top:
fileID = fopen('Top.txt','w');
fprintf(fileID,'%d\n',nodeset_top);
fclose(fileID);

% Save Bottom:
fileID = fopen('Bottom.txt','w');
fprintf(fileID,'%d\n',nodeset_bottom);
fclose(fileID);
end
