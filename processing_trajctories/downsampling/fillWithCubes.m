function [numX,numY,numZ] = fillWithCubes(minVerts, maxVerts, cubeSize)
%finds the number of cubes that can be placed along each dimension to fill
%the box as closely as possible, rounds up to ensure full space is covered
x = maxVerts(1) - minVerts(1);
y = maxVerts(2) - minVerts(2);
z = maxVerts(3) - minVerts(3);
numX = ceil(x/cubeSize);
numY = ceil(y/cubeSize);
numZ = ceil(z/cubeSize);
end