function writeObjFile(fullFileName, vertices, faces)
%fullFileName should include the path, filename and extension
fid = fopen(fullFileName, 'w');
%write info as comment at the top of the file
fprintf(fid, '# number of vertices: %d\n', size(vertices,1));
if nargin == 3
    fprintf(fid, '# number of faces: %d\n', size(faces,1));
end
%write vertices to file
for i=1:size(vertices,1)
    fprintf(fid,'%s %5.6f %5.6f %5.6f\n', 'v', vertices(i,1),vertices(i,2),vertices(i,3));
end
fprintf('\n');%new line to visually separate vertices and faces

if nargin == 3
    %write faces to file
    for j=1:size(faces,1)
        fprintf(fid,'%s %d %d %d\n', 'f', faces(j,1),faces(j,2),faces(j,3));
    end
end

fclose(fid);
end