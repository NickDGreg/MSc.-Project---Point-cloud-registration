filenamesPath = 'C:\Users\User\Documents\Imperial\msc_project\';
filenamesCSV = char(strcat(filenamesPath, 'filenames.csv'));
filenamesData = importFilenames(filenamesCSV, 2);

%used in downsampling the catheter path data or reading in the existing
%downsampled data
gridSize = 2;

% the ground truth cath data and the aorta are rotated and translated not
% the cath data beginning at the origin. This is so training data of the
% centering positions cannot be learnt simply be all aortas and
% trajectories being in the same location and orientation. It must be
% learnt from the characteristics of the catheter data

%parameters to determine the type of data created
save_data = true;
saveToOneFile = true;
inclCathStats = true;
sparse = true;
inclSimStats = false; %returns ~500 rows per trajectory so all aortas, 
                     %5 samples from each --> 32,500 (can alter difference
                     %vars %in get SimulatedTrajStats for more or less data)
oneFilename = '\Aorta Mapping Data\combined_data\cath1_sparse_GT_k3.csv';
oneFolderTrajs = '\Aorta Mapping Data\combined_data\trajectories\maxRot_sparse_';
% saveFilename = 'simData_largeTransAndRot.csv';
% trajFilename = 'trajs_simData_largeTransAndRot.csv';
saveFilename = 'simData_GT.csv';
trajFilename = 'trajs_cath1.csv';
%number of gaussians to split the aorta into (3 or 4 is good)
k = 3;
applyRandRot = false;
applyMaxRot = true;
maxRot = 30;
applyRandTranslation = false;
minTranslation = 10;
maxTranslation = 30;
%number of different data sets to create, needs to be at least one
%rotations only applied of appyRandRots is true
numRotatedSamples = 1; 

numStats = 52; %max size matrix needs to be to store all possible stats 
inclAortaMeans = true;
inclAortaEigVecs = true;
inclAortaEigVals = false;
inclAortaCov = true;

%create train data 
for row=1:13
    fprintf('Row: %i\n',row)
    if row==7
        continue
    end
    procedureFolder = char(strcat(filenamesData(row, 1)));
    procedurePath = char(strcat(filenamesPath,'Aorta Mapping Data\',procedureFolder,'\'));
	aorta_centered_Filename = char(strcat(procedurePath, filenamesData(row, 8)));
    
    saveFilepath = char(strcat(procedurePath,'\datasets\',saveFilename));
    trajFilepath = char(strcat(procedurePath,'\datasets\',trajFilename));  

    disp('reading in aorta file');
    [ aortaverts_orig, aortafaces_orig ] = read_vertices_and_faces_from_obj_file( aorta_centered_Filename );
    disp('aorta file read')
    
    disp('Get downsampled catheter trajectory')
    ds_GT_cath_filename = char(strcat(procedurePath, 'tipInfo_centered_DSgridsize2', num2str(gridSize), '.csv'));
%     ds_GT_cath_filename = char(strcat(procedurePath, 'cathMoreSim.csv'));
    try
        if sparse
            ds_GT_cath_filename = char(strcat(procedurePath, 'tipInfo_centered_sparse.csv'));
            ds_GT_cathInfo_orig = csvread(ds_GT_cath_filename);
        else
            ds_GT_cath_filename = char(strcat(procedurePath, 'tipInfo_centered_DSgridsize2', num2str(gridSize), '.csv'));
            ds_GT_cathInfo_orig = csvread(ds_GT_cath_filename);
        end
        disp('downsampled trajectory retrieved')
    catch
        fprintf(['The downsampled trajectory of that grid size does not exist.\n'...
             'So read in orignal, downsample and save downsampled version.\n'...
             'Takes a few minutes...'])

        disp('reading in catheter data');
        GT_cathInfoFilename = char(strcat(procedurePath, 'tipInfo_centered.csv'));
        GT_cathInfo = csvread(GT_cathInfoFilename);
%         cathInfoFilename = char(strcat(procedurePath, 'tipInfo_origin_start.csv'));
%         cathInfo = csvread(cathInfoFilename);
        disp('catheter data read');

        [ds_GT_cathInfo_orig, ordervisited, ~] = cubeDownsampling(GT_cathInfo, gridSize); %order visited will likely be useful when trying to plot trajectories
        csvwrite(ds_GT_cath_filename, ds_GT_cathInfo);
        disp('Downsampled and saved')
    end
    
    %get the downsampled un registered data by finding translation between
    %origin start path and ground truth, then moving the downsampled ground
    %truth by that transformation
    ds_origin_start_orig = ds_GT_cathInfo_orig;
    ds_origin_start_orig(:,1:3) = ds_origin_start_orig(:,1:3) - ds_GT_cathInfo_orig(1,1:3);
    
    disp('reading in centerline and cross section data');
    centerline_filename = char(strcat(procedurePath, 'ordered_centerline.csv'));
    centerline = csvread(centerline_filename);
    numCLpoints = size(centerline,1);
    crossSections_filename = char(strcat(procedurePath, 'crossSections.csv'));
    crossSections = csvread(crossSections_filename);    
    crossSectionsIndices_filename = char(strcat(procedurePath, 'crossSection_endIndices.csv'));
    crossSection_endIndices = csvread(crossSectionsIndices_filename);
    largestCrossSection = crossSection_endIndices(1);
    crossSection_endIndices = crossSection_endIndices(2:end);
    crossSections = reshape(crossSections,largestCrossSection,3,numCLpoints); 
    disp('centerline and cross section data read');
    
    %fit gaussian mixture model to aorta to describe shape in few
    %parameters
    [gmfit, clusterX] = getAortaDescriptors(aortaverts_orig, centerline, crossSections,...
                                            crossSection_endIndices, k);
  
    for rots=1:numRotatedSamples
        fprintf('Sample number: %d\n',rots)
        plusMinus = [1 -1];
        ds_GT_cathInfo = ds_GT_cathInfo_orig;
        aortaverts = aortaverts_orig;
        if applyRandTranslation
            %add translation to the gmfit means
            translation = randi([minTranslation maxTranslation],1,3) .* [plusMinus(randi([1 2],1)) plusMinus(randi([1 2],1)) plusMinus(randi([1 2],1))];
%             translation = [50 50 50];
            aortaMeans = gmfit.mu + translation;
            aortaverts = aortaverts + translation;
%             centerline = centerline + translation;
            ds_GT_cathInfo(:,1:3) = ds_GT_cathInfo(:,1:3) + translation;
        else 
            translation = [0 0 0];
            aortaMeans = gmfit.mu;
        end
        if applyRandRot
            %rotate gmfit data cov, extract the
%             rng(5);
            xAngle =rand * maxRot * plusMinus(randi([1 2],1));
%             rng(3);
            yAngle =rand * maxRot * plusMinus(randi([1 2],1));
%             rng(1);
            zAngle =rand * maxRot * plusMinus(randi([1 2],1));
            rotMat = createRotationMatrix(xAngle,yAngle,zAngle);
        elseif applyMaxRot
            rotMat = createRotationMatrix(maxRot,maxRot,maxRot);
        else
            rotMat = eye(3);
        end
        aortaverts = rotatePositions(aortaverts, rotMat);
%        centerline = rotatePositions(centerline, rotMat);
        aortaMeans = rotatePositions(aortaMeans, rotMat);
        ds_GT_cathInfo(:,1:3) = rotatePositions(ds_GT_cathInfo(:,1:3), rotMat);
        
        %get the cathetory trajectory data
        if inclCathStats
            disp('getting the catheter trajectory stats...')
            [ cathStats, statNames, statEndIdxs, trajectoriesMat ] = getCathTrajStats( ds_origin_start_orig, ds_GT_cathInfo, numStats, row, rots );
            cathStats(:,1) = row;
            cathStats(:,2) = 1; %signal cath data or simulated data
            disp('done')
        end

        %get the simulated trajectory data
        if inclSimStats
            disp('getting the simulated trajectory stats...')
            [ simulatedStats, simulatedTrajs, statNames_, statEndIdxs_ ] = getSimulatedTrajStats( procedurePath, translation, rotMat, numStats, row, rots );
            simulatedStats(:,1) = row;
            simulatedStats(:,2) = 0; %signal cath data or simulated data
            if inclCathStats && row==1 
                cathStats_old = cathStats;
                cathStats = [cathStats(1,:);
                             simulatedStats(1,:);
                             cathStats(2:end,:);
                             simulatedStats(2:end,:)];
            elseif inclCathStats
                cathStats_old = cathStats;
                cathStats = [cathStats;
                             simulatedStats];
            else
                statNames = statNames_;
                statEndIdxs = statEndIdxs_;
                cathStats = simulatedStats;
            end
            disp('done')
        end
    
        aortaDescriptors = ones(1,18*k);
        aortaDescriptorIdxs = ones(1,4*k);
        aortaDesc_idx = 0;
        numAortaDesc = 0;
        if inclAortaMeans
            aortaDescriptors(aortaDesc_idx+1:aortaDesc_idx+k*3) = reshape(aortaMeans',1,k*3);
            aortaDescriptorIdxs(numAortaDesc+1:numAortaDesc+k) = (1:k)*k;
            kNames = textscan(num2str(1:k),'%s');
            meanCell = cell(1,k);
            meanCell(:) = {'mean'};
            aortaDescriptorNames(numAortaDesc+1:numAortaDesc+k) = cellfun(@(x,y) [x '_' y],meanCell,kNames{1}','un',0);
            aortaDesc_idx = aortaDesc_idx+k*3;
            numAortaDesc = numAortaDesc + k;
        end
        
        %add the aorta descriptors if they're selected
        if inclAortaEigVecs || inclAortaEigVals || inclAortaCov
            for covMat=1:k
                rotCov = rotMat * gmfit.Sigma(:,:,covMat) * rotMat';
                [D,V] = eig(rotCov);
                if inclAortaEigVecs
                    aortaDescriptors(aortaDesc_idx+1:aortaDesc_idx+9) = reshape(D,1,9); 
                    aortaDescriptorIdxs(numAortaDesc+1) = aortaDesc_idx+9;
                    aortaDesc_idx = aortaDesc_idx+9;
                    aortaDescriptorNames{numAortaDesc+1} = strcat('aortaEigVec_', num2str(covMat));
                    numAortaDesc = numAortaDesc + 1;
                end
                
                if inclAortaEigVals
                    aortaDescriptors(aortaDesc_idx+1:aortaDesc_idx+3) = diag(V); 
                    aortaDescriptorIdxs(numAortaDesc+1) = aortaDesc_idx+3;
                    aortaDesc_idx = aortaDesc_idx+3;
                    aortaDescriptorNames{numAortaDesc+1} = strcat('aortaEigVal_', num2str(covMat));
                    numAortaDesc = numAortaDesc + 1;
                end
                
                if inclAortaCov
                    aortaDescriptors(aortaDesc_idx+1:aortaDesc_idx+3) = rotCov(logical(tril(ones(3),-1))); %lower triangular of cov mat
                    aortaDescriptorIdxs(numAortaDesc+1) = aortaDesc_idx+3;
                    aortaDesc_idx = aortaDesc_idx+3;
                    aortaDescriptorNames{numAortaDesc+1} = strcat('aortaCov_', num2str(covMat));
                    numAortaDesc = numAortaDesc + 1;
                end
            end
        end
        
        if numAortaDesc > 0
            aortaDescriptors = aortaDescriptors(1:aortaDesc_idx);
            %insert aorta descrpitors into the cath stats
            statRow = size(cathStats,1);
            cathStats = [cathStats(:,1:end-3) repmat(aortaDescriptors,statRow,1)  cathStats(:,end-2:end)];
            
            if row==1 && rots==1
                aortaDescriptorIdxs = aortaDescriptorIdxs(1:numAortaDesc);
                aortaDescriptorNames = aortaDescriptorNames(1:numAortaDesc);

                %shift end idx to account for inserted descriptors
                statEndIdxs(end) = statEndIdxs(end) + aortaDesc_idx; 
                %shift aorta desc idxs as they are inserted after the cath stats
                aortaDescriptorIdxs = aortaDescriptorIdxs + statEndIdxs(end-1); 
                statEndIdxs = [statEndIdxs(1:end-1) aortaDescriptorIdxs statEndIdxs(end)];
                statNames = [statNames(1:end-1) aortaDescriptorNames statNames(end)];
            end
        end
        
        %alter this because theres no way its as good as it could be 
        if save_data
            if saveToOneFile
                oneFilepath = char(strcat(filenamesPath,oneFilename));
                trajFilepath = char(strcat(filenamesPath,oneFolderTrajs,...
                                    strcat('traj_',int2str(row),'.csv')));
                if row==1 && rots == 1
                    filePh = fopen(oneFilepath,'w');
                    for header=1:2
                        for i=1:numStats-1 + (k*18) 
                            try
                                if header==1
                                    fprintf(filePh,'%s,',statNames{i});
                                else
                                    fprintf(filePh,'%d,',statEndIdxs(i));
                                end
                            catch
                                fprintf(filePh,',');
                            end
                        end
                        fprintf(filePh,'\n');
                    end
                    fclose(filePh);
                    dlmwrite(oneFilepath,cathStats,'-append','delimiter', ',')
                    csvwrite(trajFilepath, trajectoriesMat)
                else
                    dlmwrite(oneFilepath,cathStats,'-append','delimiter', ',')
                    csvwrite(trajFilepath, trajectoriesMat)
                end
            else
                if row==1 && rots == 1
                    filePh = fopen(saveFilepath,'w');
                    for header=1:2
                        for i=1:numStats-1 + (k*18) 
                            try
                                if header==1
                                    fprintf(filePh,'%s,',statNames{i});
                                else
                                    fprintf(filePh,'%d,',statEndIdxs(i));
                                end
                            catch
                                fprintf(filePh,',');
                            end
                        end
                        fprintf(filePh,'\n');
                    end
                    fclose(filePh);
%                     dlmwrite(saveFilepath,cathStats,'-append','delimiter', ',')
                    dlmwrite(trajFilepath,trajectoriesMat,'-append','delimiter', ',')
                else
%                     dlmwrite(saveFilepath,cathStats,'-append','delimiter', ',')
                    dlmwrite(trajFilepath,trajectoriesMat,'-append','delimiter', ',')
                end
            end
        end
    end
end
% % Project the scan, catheters so that it aligns with the reg initialising centers
% % found in the stats mat which uses PCA, otherwise plot makes no sense
% % 
% [ aortavertsOrig, aortafacesOrig ] = read_vertices_and_faces_from_obj_file( aorta_centered_Filename );
% 
% figure, hold on, view(3), grid on;
% % patch('Faces',aortafacesOrig,'Vertices',aortaverts, 'EdgeColor', 'blue', 'EdgeAlpha', 0.05, 'FaceColor', 'none');
% scatter3(aortaverts(:,1),aortaverts(:,2),aortaverts(:,3),2,'blue','filled');
% scatter3(cathStats(:,end-2),cathStats(:,end-1),cathStats(:,end),4,'red','filled');
% scatter3(ds_GT_cathInfo(:,1),ds_GT_cathInfo(:,2),ds_GT_cathInfo(:,3),1,'green','filled');
% scatter3(aortaMeans(:,1),aortaMeans(:,2),aortaMeans(:,3));
% % aortaMeans_ = reshape(aortaMeans,k,3)';
% % aortaMeans_ = [aortaMeans_(:,1:3);
% %                aortaMeans_(:,4)'];
% % aortaMeans_ = reshape(aortaMeans,3,k)';
% % aortaMeans_ = [repmat(aortaMeans_(1,:),3,1);
% %                repmat(aortaMeans_(2,:),3,1);
% %                repmat(aortaMeans_(3,:),3,1)];
% %            
% % aortaDescriptors_ = reshape(aortaDescriptors,15,k)';
% % aortaEigVecs = aortaDescriptors_(:,1:9);
% % aortaEigVecs_ = reshape(aortaEigVecs',3,12)';
% % aortaEigVs = aortaDescriptors_(:,10:12);
% % quiver3(aortaMeans_(:,1),aortaMeans_(:,2),aortaMeans_(:,3),aortaEigVecs_(:,1),aortaEigVecs_(:,2),aortaEigVecs_(:,3))
% % scatter3(aortaMeans_(:,1),aortaMeans_(:,2),aortaMeans_(:,3));
% % 
% % patch('Faces',aortafacesOrig,'Vertices',aortavertsOrig, 'EdgeColor', 'red', 'EdgeAlpha', 0.05, 'FaceColor', 'none');
% % aortaMeans_orig = gmfit.mu;
% % aortaDescriptors_orig = ones(1,15*k);
% % for covMat=1:k
% %     start = (15*(covMat-1))+1;
% %     rotCov = rotMat * gmfit.Sigma(:,:,covMat) * rotMat';
% %     [D,V] = eig(rotCov);
% %     aortaDescriptors_orig(start:start+8) = reshape(D,1,9); %eigenvectors
% %     aortaDescriptors_orig(start+9:start+11) = diag(V); %eigenvalues
% %     aortaDescriptors_orig(start+12:start+14) = rotCov(logical(tril(ones(3),-1))); %lower triangular of cov mat
% % end
% % 
% % aortaMeans_orig = [repmat(aortaMeans_orig(1,:),3,1);
% %                repmat(aortaMeans_orig(2,:),3,1);
% %                repmat(aortaMeans_orig(3,:),3,1)];
% %            
% % aortaDescriptors_orig = reshape(aortaDescriptors_orig,15,k)';
% % aortaEigVecs_orig = aortaDescriptors_orig(:,1:9);
% % aortaEigVecs_orig = reshape(aortaEigVecs_orig',3,12)';
% % quiver3(aortaMeans_orig(:,1),aortaMeans_orig(:,2),aortaMeans_orig(:,3),aortaEigVecs_orig(:,1),aortaEigVecs_orig(:,2),aortaEigVecs_orig(:,3))
% % scatter3(aortaMeans_orig(:,1),aortaMeans_orig(:,2),aortaMeans_orig(:,3),1,'blue','filled');
% % scatter3(cathStats(:,end-2),cathStats(:,end-1),cathStats(:,end),3,'red','filled');
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold off
% 
% fclose('all')