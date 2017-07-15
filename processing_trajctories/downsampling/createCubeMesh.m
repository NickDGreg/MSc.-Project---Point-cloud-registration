function [vertices, faces] = createCubeMesh(vmin, vmax)
%vmin is [x y z] minimum values of the vertices, vmax max values
vertices = [vmax(1) vmin(2) vmin(3); 
            vmax(1) vmax(2) vmin(3);   
            vmin(1) vmax(2) vmin(3); 
            vmin(1) vmax(2) vmax(3); 
            vmin(1) vmin(2) vmax(3); 
            vmax(1) vmin(2) vmax(3); 
            vmin; vmax ];
        
faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];