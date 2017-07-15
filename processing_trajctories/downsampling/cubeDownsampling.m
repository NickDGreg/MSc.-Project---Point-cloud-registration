function [ds_cathInfo, orderVisited, cubes, cubes_plot] = cubeDownsampling(cathInfo, cubeSize)
%input (points, cubessize, ture/false) true or empty indicates catheter
%data
%downsamples trajectory by creating cuboid around first point, then averaging
%all points within cube, if point is outside, new cube is appended to the
%existing cubes, by using a data structure that corresponds
%to the distances can directly index the cube rather than searching through
%a list. SLOW though, ~45-1min for 10,000 points, but real time, maybe in a
%faster language it may be ok

numPoints = size(cathInfo,1);
numFields = size(cathInfo,2);
pointAssignment = zeros(numPoints,3); %not currently used, but maybe useful

vmin = min(cathInfo(:,1:3));
vmax = max(cathInfo(:,1:3));
[numX,numY,numZ] = fillWithCubes(vmin, vmax, 1);
%create 3D array of uninitialised cubes, dimensions fit the space as closely 
%as possible using user defined cube size
cubes = repmat( struct('minVerts',[0 0 0],'maxVerts',[0 0 0],'visited',0,'average',zeros(1,numFields),'merged',0),numX,numY,numZ);
cubesVisited = zeros(numPoints, 3); %list to record cubes that are moved through (initialised)
newCube = 0; %index for how many cubes are instantiated
orderVisited = zeros(numPoints,1); %list to track how the catheter moves through the cubes
visited = 0; %index to track the order the cubes are visited in 

%set first cube with first point at center
p = cathInfo(1,1:3);
iCube = ceil( (p - vmin) / cubeSize ); %get index of point using point loc and min x,y,z of all points
iCube(~iCube) = 1; %ensure zeros become 1s (0 if point is the minimum x,y,or z)
cubes(iCube(1),iCube(2),iCube(3)).minVerts = vmin + (iCube-1)*cubeSize;
cubes(iCube(1),iCube(2),iCube(3)).maxVerts = vmin + (iCube)*cubeSize;
cubes(iCube(1),iCube(2),iCube(3)).visited = 1;
cubes(iCube(1),iCube(2),iCube(3)).average = cubes(iCube(1),iCube(2),iCube(3)).average + cathInfo(1,:);
cubes(iCube(1),iCube(2),iCube(3)).merged = 1;
visited = visited + 1;
newCube = newCube + 1;
orderVisited(visited) = newCube;
cubesVisited(newCube,:) = iCube;

%uncomment for real time plot
% [verts, faces] = createCubeMesh(cubes(iCube(1),iCube(2),iCube(3)).minVerts, cubes(iCube(1),iCube(2),iCube(3)).maxVerts); 
% figure, hold on, view(3), grid on
% plot3(p(1),p(2),p(3), 'r.', 'MarkerSize', 5);
% patch('Vertices',verts,'Faces',faces,'FaceColor','green','FaceAlpha',0.2);
% % patch('Vertices',mainCubeVerts,'Faces',mainCubeFaces,'FaceColor','None');

tic
for i=2:numPoints
    p = cathInfo(i,1:3);
%     plot3(p(1),p(2),p(3), 'r.', 'MarkerSize', 5);   %uncomment for real time plot
%     drawnow    %uncomment for real time plot
%     pause(0.000001) %uncomment for real time plot
    
    %check if point in same cube as previous point, thats most likely with catheter data 
    if pointInCube(cubes(iCube(1),iCube(2),iCube(3)).minVerts, cubes(iCube(1),iCube(2),iCube(3)).maxVerts, p)
        cubes(iCube(1),iCube(2),iCube(3)).average = cubes(iCube(1),iCube(2),iCube(3)).average + cathInfo(i,:);
        cubes(iCube(1),iCube(2),iCube(3)).merged = cubes(iCube(1),iCube(2),iCube(3)).merged + 1;
        pointAssignment(i,:) = iCube;
    else %locate cube
        iCube = ceil( (p - vmin) / cubeSize );
        iCube(~iCube) = 1; %ensure zero indices become 1s
        if cubes(iCube(1),iCube(2),iCube(3)).visited
            %safety check for failures
            if pointInCube(cubes(iCube(1),iCube(2),iCube(3)).minVerts, cubes(iCube(1),iCube(2),iCube(3)).maxVerts, p)
                cubes(iCube(1),iCube(2),iCube(3)).average = cubes(iCube(1),iCube(2),iCube(3)).average + cathInfo(i,:);
                cubes(iCube(1),iCube(2),iCube(3)).merged = cubes(iCube(1),iCube(2),iCube(3)).merged + 1;
                pointAssignment(i,:) = iCube;
                locationArray = ismember(cubesVisited, iCube, 'rows');
                loc = find(locationArray, 1);
                visited = visited + 1;
                orderVisited(visited) = loc;
            else
                disp('Point isnt in the cube that its index suggests?!')
                p
                iCube 
                minVerts = cubes(iCube(1),iCube(2),iCube(3)).minVerts
                maxVerts = cubes(iCube(1),iCube(2),iCube(3)).maxVerts
            end
        else
            %new cube so initialise cube 
            cubes(iCube(1),iCube(2),iCube(3)).minVerts = vmin + (iCube-1)*cubeSize;
            cubes(iCube(1),iCube(2),iCube(3)).maxVerts = vmin + (iCube)*cubeSize;
            cubes(iCube(1),iCube(2),iCube(3)).visited = 1;
            cubes(iCube(1),iCube(2),iCube(3)).average = cubes(iCube(1),iCube(2),iCube(3)).average + cathInfo(i,:);
            cubes(iCube(1),iCube(2),iCube(3)).merged = 1;
            pointAssignment(i,:) = iCube;
            visited = visited + 1;
            newCube = newCube + 1;
            orderVisited(visited) = newCube;
            cubesVisited(newCube,:) = iCube;
          
            %uncomment for real time plot
%             [verts, faces] = createCubeMesh(cubes(iCube(1),iCube(2),iCube(3)).minVerts, cubes(iCube(1),iCube(2),iCube(3)).maxVerts);
%             patch('Vertices',verts,'Faces',faces,'FaceColor',rand(1,3),'FaceAlpha',0.2);
        end
    end
end
toc
cubesVisited = cubesVisited(1:newCube,:);
orderVisited = orderVisited(1:visited);

tic
disp('collecting downsampled data')
ds_cathInfo = zeros(newCube, numFields);
cubes_plot = zeros(newCube, 6); %save to plot how averages were taken

for i=1:newCube
    iCube = cubesVisited(i,:);
    cubes(iCube(1),iCube(2),iCube(3)).average = cubes(iCube(1),iCube(2),iCube(3)).average / cubes(iCube(1),iCube(2),iCube(3)).merged; 
    ds_cathInfo(i, :) = cubes(iCube(1),iCube(2),iCube(3)).average;
    cubes_plot(i,:) = [cubes(iCube(1),iCube(2),iCube(3)).minVerts cubes(iCube(1),iCube(2),iCube(3)).maxVerts]; %save to plot how averages were taken
end



toc
end