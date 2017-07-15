datafolder = 'C:\Users\User\Documents\Imperial\msc_project\Aorta Mapping Data\';
filenamesPath = 'C:\Users\User\Documents\Imperial\msc_project\';
filenamesCSV = char(strcat(filenamesPath, 'filenames.csv'));

%requires manual observation to determine if the heuristic method to remove
%the datapoints from heart chamber exploration is successful, some
%procedures enter heart chambers explore, re-enter aorta and then explore
%again. Hence may need to remove heart chamber exploration twice.

filenames = importFilenames(filenamesCSV, 2);
numRows = size(filenames,1);
origin = [0 0 0];

%only select one of these at a time, make sure the correct cath points are
%selected for complete
firstRun = 0;
rerun = 0;
complete = 1;
%can be added with rerun
moreThanTwoRuns = 0;
%add with complete once checked correct alignment plot
saveCath = 0;


tipInfo_aligned_file = char(strcat(path, 'tipInfo_Aligned.obj')); 
[ aortaTipposII_align, ~ ] = read_vertices_and_faces_from_obj_file( tipInfo_aligned_file );

% If manual registration was required to correct original for ground truth
% alter the orientation by the rotation matrix
% lastIndexII = 17402;
% aortaTipposII = aortaTipposII_align(1:lastIndexII,:);
% % aortaOrientII = aortaOrient(1:lastIndexII,:);
% aortaOrientII = tiporient(1:lastIndexII,:);
% transformationMatrix = [1.00    -0.05	-0.01;	
%                         0.05	1.00    0.00;	
%                         0.01	0.00   1.00];	
% aortaOrientII = aortaOrientII*transformationMatrix';

for row=13:13
    if rerun && firstRun
        error('saftey check, cannot be the first run and a rerun at the same time')
    end
    
    if firstRun
        procedureFolder = char(strcat(filenames(row, 1)));
        path = char(strcat(datafolder,procedureFolder,'\'));
        niobefilename = char(strcat(path, filenames(row, 2)));
        aortafilename = char(strcat(path, filenames(row, 5))); %original aorta file
        aortafilename_centered = char(strcat(path, filenames(row, 8))); %centered aorta file
        
        disp('reading in aorta file');
        [ aortaverts, aortafaces ] = read_vertices_and_faces_from_obj_file( aortafilename );    
        disp('reading in niobe file');
        [ tippos_orig, tiporient_orig, cartotime, contact] = readnewfilteredniobefile( niobefilename );
        %automatically cut most off which will be in the heart chambers
        %speeds up the process
        tippos = tippos_orig(1:200000,:);
        tiporient = tiporient_orig(1:200000,:);

        disp('calculating points in aorta');
        [in, pointsin] = pointsinaorta( aortaverts, aortafaces, tippos);
        disp('found points in aorta')

        cartotime_in = cartotime(in);

    %     pointsin - human registered aorta from only the first path up the aorta
    %     aortaInOut bool if point in first path is in/out the human reg aorta
    %     aortaTippos all positions catheter moved through during first path
        disp('get aorta path, exclude heart data')
        [pointsin, aortaInOut, aortaTippos, aortaOrient, lastIndex] = firstAortaPathPoints(in, cartotime_in, cartotime, pointsin, tippos, tiporient);
        shortened_cathInfo = [aortaTippos aortaOrient contact(1:lastIndex)];

    %     test plot of aorta to check if need to re-run heart chamber exploration
        colouration = 1:size(shortened_cathInfo,1);
        figure, hold on, view(3);
        patch('Faces',aortafaces,'Vertices',aortaverts, 'EdgeColor', 'blue', 'EdgeAlpha', 0.03, 'FaceColor', 'none');
        scatter3(shortened_cathInfo(:,1), shortened_cathInfo(:,2), shortened_cathInfo(:,3), 2, colouration', 'filled'); 
        colormap(jet)
        colorbar()
    %     scatter3(pointsin(:,1), pointsin(:,2), pointsin(:,3), 1, 'green', 'filled'); 
    %     scatter3(aortaTippos(~aortaInOut,1), aortaTippos(~aortaInOut,2), aortaTippos(~aortaInOut,3), 1, 'blue', 'filled');
    %     scatter3(tippos(:,1), tippos(:,2), tippos(:,3), 1, 'blue', 'filled');
        legend('aorta', 'in', 'out')
        hold off

        lastIndex
    end

%     re run if it still captures a lot of heart mapping
    if rerun
        disp('re run aorta path check')
        if moreThanTwoRuns
            cartotime_inII = cartotime(aortaInOutII);
            [pointsinII, aortaInOutII, aortaTipposII, aortaOrientII, lastIndexII] = firstAortaPathPoints(aortaInOutII, cartotime_inII, cartotime, pointsinII, aortaTipposII, aortaOrientII);
        else
            cartotime_inII = cartotime(aortaInOut);
            [pointsinII, aortaInOutII, aortaTipposII, aortaOrientII, lastIndexII] = firstAortaPathPoints(aortaInOut, cartotime_inII, cartotime, pointsin, aortaTippos, aortaOrient);
        end
        lastIndexII
        colouration = (1:size(aortaTipposII,1))';
    %   %test plot of aorta to check if need to re-run heart chamber exploration
        figure, hold on, view(3);
        patch('Faces',aortafaces,'Vertices',aortaverts, 'EdgeColor', 'red', 'EdgeAlpha', 0.03, 'FaceColor', 'none');
        scatter3(aortaTipposII(:,1), aortaTipposII(:,2), aortaTipposII(:,3), 1, colouration, 'filled');
        colormap(jet)
        colorbar()
        legend('aorta', 'in', 'out')
        hold off
    end
end


if complete
    %align the catheter points with the centered aorta
    [aortaverts_centered, aortaDiff] = alignCenters(origin, aortaverts);
    cathAligned = aortaTipposII + aortaDiff;
    
    %find start Index
    current = cathAligned(1,:);
    for i=2:size(cathAligned,1)
        next = cathAligned(i,:);
        if next == current
            current = next;
        else
            startIndex = i;
            break;
        end
    end
    
    colouration = (1:size(aortaTipposII,1))';
%   %test plot to check if points are correctly aligned
    figure, hold on, view(3);
    patch('Faces',aortafaces,'Vertices',aortaverts_centered, 'EdgeColor', 'red', 'EdgeAlpha', 0.03, 'FaceColor', 'none');
% 	scatter3(cathAligned(:,1), cathAligned(:,2), cathAligned(:,3), 1, colouration, 'filled');
    quiver3(cathAligned(:,1),cathAligned(:,2),cathAligned(:,3),aortaOrientII(:,1),aortaOrientII(:,2),aortaOrientII(:,3))
    colormap(jet)
    colorbar()
    legend('aorta', 'in', 'out')
    hold off
    
    finalCathInfo = [cathAligned aortaOrientII contact(1:lastIndexII)];
    %save final file
    cathInfo_filename = char(strcat(path, 'tipInfo_centered.csv'));
    finalCathInfo = [cathAligned aortaOrientII contact(1:lastIndexII)];
    if saveCath
        csvwrite(cathInfo_filename, finalCathInfo);
    end
end

% %If misalinged, save obj file, realign in meshlab then save
% colouration = (1:size(aortaTipposII,1))';
% % %   %test plot of aorta to check if need to re-run heart chamber exploration
% figure, hold on, view(3);
% patch('Faces',aortafaces,'Vertices',aortaverts, 'EdgeColor', 'red', 'EdgeAlpha', 0.03, 'FaceColor', 'none');
% scatter3(aortaTipposII(:,1), aortaTipposII(:,2), aortaTipposII(:,3), 1, colouration, 'filled');
% colormap(jet)
% colorbar()
% legend('aorta', 'in', 'out')
% hold off
% cathInfo_misName = char(strcat(path, 'tipInfo_misAligned.obj'));
% cathInfo_misaligned = [aortaTipposII aortaOrientII contact(1:lastIndexII)];
% writeObjFile(cathInfo_misName, aortaTippos);