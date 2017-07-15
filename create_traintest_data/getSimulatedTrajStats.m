function [ simulatedStats, simulatedTrajs, statNames, statEndIdxs ] = getSimulatedTrajStats( procedurePath, translation, rotMat, numStats, row, rots )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%alter the density at the start or end of the path to simulate more points
%being collected in one region
diffDensities = 1;

%all as percentages of total num points
earliestStart = 0.2;
minSampleSize = 0.03;
latestStart = 0.6; %no subset can start later than this percentage of the points
sampleSteps = 0.03; %difference between starting positions of samples of the same size
sampleSizeDiff = 0.03; % n+1 set of samples contain ~3% of total cath points more than n set
peturbation = 0.00; %random perturbation of start and end of the samples 
numRandRows = ceil(1/sampleSteps * 1/sampleSizeDiff);


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

disp('create the orientations')
orientations = zeros(numCLpoints, 3);
for CLpoint=2:numCLpoints
    CLvec = centerline(CLpoint,:)-centerline(CLpoint-1,:);
    %calculate the direction cosine
    rotations = [CLvec(1)/norm(CLvec) CLvec(2)/norm(CLvec) CLvec(3)/norm(CLvec)];
    orientations(CLpoint,:) = rotations; 
end   
orientations(1,:) = orientations(2,:);
disp('done')

%get the stats matrices for the subsets of the simulated data
trajCol = 1;
trajectoriesMat = zeros(300, 8000);
longestTraj = 0;
simulatedStats = zeros(numRandRows,numStats); %last column will be blank, only filled in not random case, see catheterPathStats.m
statRow = 0;
lastStartIndex = ceil(latestStart*numCLpoints);

if row==1 && rots==1
    simulatedStats(1,1:5) = [latestStart minSampleSize sampleSteps sampleSizeDiff peturbation];
    statRow = statRow + 1;
end

disp('get the simulated stats')
for gap=minSampleSize:sampleSizeDiff:1
    peturb = ceil(gap*numCLpoints*peturbation);
    for p=earliestStart:sampleSteps:1-gap  
        start = ceil(p*numCLpoints);
        if p+gap < numCLpoints && start<=lastStartIndex
            start = ceil(p*numCLpoints);
            sidx = max( [start+randi([-peturb,peturb],1) 1] );
            eidx = min( [start+ceil(gap*numCLpoints)+randi([-peturb,peturb],1) numCLpoints]);
            midpoint = ceil((sidx+eidx)/2);

            statRow = statRow + 1;
            orientationSidx = repmat(orientations(sidx,:),crossSection_endIndices(sidx),1);
            orientationMidx = repmat(orientations(midpoint,:),crossSection_endIndices(midpoint),1);
            orientationEidx = repmat(orientations(eidx,:),crossSection_endIndices(eidx),1);
            simulatedData = [centerline(sidx:eidx,:) orientations(sidx:eidx,:); 
                             crossSections(1:crossSection_endIndices(sidx),:,sidx) orientationSidx;
                             crossSections(1:crossSection_endIndices(midpoint),:,midpoint) orientationMidx;
                             crossSections(1:crossSection_endIndices(eidx),:,eidx) orientationEidx];
            simulatedStats(statRow,3:4) = [sidx eidx];
            simulatedStats(statRow,3:4) = [sidx eidx];
            
            %translate simulated data to origin and transform
            %groundtruth
            simulatedData_origin = simulatedData - simulatedData(1,:);
            simulatedData(:,1:3) = simulatedData(:,1:3) + translation;
            simulatedData(:,1:3) = rotatePositions(simulatedData(:,1:3),rotMat);
            simulatedData(:,4:6) = rotatePositions(simulatedData(:,4:6), rotMat);

            [simStatsRow,~,statEndIdxs] = getStatsMatRow(simulatedData_origin, simulatedData);
            simulatedStats(statRow,5:statEndIdxs(end)+4) = simStatsRow;
            
            %save the trajectory
            trajEidx = size(simulatedData,1);
            trajectoriesMat(2,trajCol) = trajEidx;
            trajectoriesMat(3:trajEidx+2,trajCol:trajCol+2) = simulatedData(:,1:3);
            trajCol = trajCol + 3;
            if trajEidx > longestTraj
               longestTraj = trajEidx;
            end

            %alter the densities if the set of catheter points is
            %spread out along the aorta, else there's no point as
            %all the rings are close together
            dataSpread = (eidx-sidx)/eidx;
            if diffDensities && dataSpread > 0.3
%                 if diffDensities && (p+gap)*numCLpoints > 0.8
                sections = [sidx midpoint eidx];
                for density=6:4:10
                    for sectionIdx=1:2:3
                        statRow = statRow + 1;
                        section = sections(sectionIdx);                                
                        otherSections = sections(sections ~= section);
                        %crossSection_endIndices is index for the
                        %last point in cross section which is not
                        %zero (cross section size had to be
                        %predefined to deal with varying sizes)
                        repeat = [crossSections(1:crossSection_endIndices(section),:,section) repmat(orientations(section,:),crossSection_endIndices(section),1)];
                        repeatedSection = repmat(repeat,density,1);
                        simulatedDataRepeat = [centerline(sidx:eidx,:) orientations(sidx:eidx,:); 
                                         repeatedSection;
                                         crossSections(1:crossSection_endIndices(otherSections(1)),:,otherSections(1)) repmat(orientations(otherSections(1),:),crossSection_endIndices(otherSections(1)),1);
                                         crossSections(1:crossSection_endIndices(otherSections(2)),:,otherSections(2)) repmat(orientations(otherSections(2),:),crossSection_endIndices(otherSections(2)),1)];
                        simulatedStats(statRow,3:4) = [sidx eidx];

                        %transform data to origin and apply
                        %transformation for the ground truth
                        simulatedDataRepeat_origin = simulatedDataRepeat - simulatedDataRepeat(1,:);
                        simulatedDataRepeat(:,1:3) = simulatedDataRepeat(:,1:3) + translation;
                        simulatedDataRepeat(:,1:3) = rotatePositions(simulatedDataRepeat(:,1:3),rotMat);
                        simulatedDataRepeat(:,4:6) = rotatePositions(simulatedDataRepeat(:,4:6), rotMat);
                        [simStatsRow,statNames,statEndIdxs] = getStatsMatRow(simulatedDataRepeat_origin, simulatedDataRepeat);
                        simulatedStats(statRow,5:statEndIdxs(end)+4) = simStatsRow;

                        %save the trajectory
                        trajEidx = size(simulatedData,1);
                        trajectoriesMat(2,trajCol) = trajEidx;
                        trajectoriesMat(3:trajEidx+2,trajCol:trajCol+2) = simulatedData(:,1:3);
                        trajCol = trajCol + 3;
                    end
                end
            end
            
        else
            break
        end
    end
end

trajCol = trajCol-3;
trajectoriesMat(1,1:7) = [trajCol, earliestStart, minSampleSize, latestStart, sampleSteps, sampleSizeDiff, peturbation];
simulatedTrajs = trajectoriesMat(1:longestTraj+2,1:trajCol+2); 
simulatedStats = simulatedStats(1:statRow,1:statEndIdxs(end)+4);
end

