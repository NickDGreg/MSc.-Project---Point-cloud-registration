function [centerlineOutput, minCLIdx, maxCLIdx] = getCenterlineStats(GT_cathSubset,centerline)
%finds the stats describing the registration centering position in terms of
%the centerline.
%centeringPosition - xyz coords of where catheter points align with aorta
%GT_cathSubset - full set of points in the subset
%Centerline - the points of the centerline of the aorta
%CLPctStep - size of the step size taken to downsample the centerline
%output: a vector containing:
            %CL point nearest to centering pos
            %Two CL points which describe the extremes of the catheter path
            %Those two points with percentages describing their indices
%minCLIdx, maxCLIdx the indices of the two extrema in the centerline

center = mean(GT_cathSubset);
nIdx = dsearchn(centerline,center);
nearest_CLPoint = centerline(nIdx,:); %nearestcenterline point to the GT centering position

maxes = max(GT_cathSubset);
mins = min(GT_cathSubset);
maxSubset = GT_cathSubset ./ ( ones(size(GT_cathSubset)) .* maxes ); %max values become 1
minSubset = GT_cathSubset ./ ( ones(size(GT_cathSubset)) .* mins ); %min values become 1
maxIdxs = find( maxSubset==1 ); %find the indices of the ones
minIdxs = find( minSubset==1 );
[Xmax,~] = ind2sub(size(GT_cathSubset), maxIdxs);
[Xmin,~] = ind2sub(size(GT_cathSubset), minIdxs);
k = dsearchn(centerline,[GT_cathSubset(Xmax,:);... %indices of the max/min values
                         GT_cathSubset(Xmin,:)]); 
minCLIdx = min(k);
maxCLIdx = max(k);
centerline(k,:);
minCLPoint = centerline(minCLIdx,:);
maxCLPoint = centerline(maxCLIdx,:);

numPoints = size(centerline,1);
p = 1/numPoints;
minIdx_as_pct = minCLIdx*p;
maxIdx_as_pct = maxCLIdx*p;

centerlineOutput = [nearest_CLPoint minCLPoint maxCLPoint minIdx_as_pct maxIdx_as_pct]; 
end