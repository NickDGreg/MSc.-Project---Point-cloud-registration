function ordered_centerline = order_centerline(centerline, threshold,maxDist)
    %Order the centerline from the tip of the descending aorta to the aortic root
    %threshold needs to be chosen based on how the sequence appears
    %different size branches require different values
    %~15 good for most, large branches (aorta 2) require 30 
    %could be faster but this is only for a small number of points 
    %and only needs to be done once 

    ordered_centerline = Inf(size(centerline)); %as checking distances, cant be zero as that could count as a nearest point
    sequence = zeros(size(centerline));
    numPoints = size(centerline,1);
    point = 1;
    idx = 1;
    seqStrtIdx = 1;
    
    %start at the tip of the descending aorta (coords with max z)
    maxes = max(centerline);
    [x, ~, ~] = ind2sub(size(centerline),find(centerline==maxes(3)));
    sequence(1,:) = centerline(x,:);
    centerline(x,:) = [];
    filledTo = 0;
    while point <= numPoints && ~isempty(centerline)
        [loc, dist] = dsearchn(centerline, sequence(idx,:));
        if dist < threshold
            idx = idx + 1;
            sequence(idx,:) = centerline(loc, :);
            centerline(loc, :) = [];
            
        elseif dist > maxDist
            disp('distance too large, must be back tracking to offshoot')
            offShootIdx = idx;
            break 
        else
%             fprintf('Branch. Distance is %d. Start index is %d \n',dist,seqStrtIdx);
            if filledTo == 0
%                 fprintf('Adding first sequence of length %d \n', idx)
                ordered_centerline(1:idx,:) = sequence(1:idx, :);
                filledTo = idx;
            else
                [offShootIdx, ~] = dsearchn(sequence,centerline(loc,:));
%                 fprintf('Adding sequence of length %d in at %d \n',idx,seqStrtIdx)
                ordered_centerline = [ordered_centerline(1:seqStrtIdx-1,:); 
                                      sequence(1:offShootIdx, :);
                                      ordered_centerline(seqStrtIdx:filledTo,:)];
                filledTo = filledTo + offShootIdx;
            end
            [seqStrtIdx, ~] = dsearchn(ordered_centerline, centerline(loc,:));
            seqStrtIdx = seqStrtIdx + 1;
            sequence = zeros(numPoints-filledTo, 3);
%             fprintf('Number of remaining points is %d \n',numPoints-filledTo);
            sequence(1,:) = centerline(loc, :);
            centerline(loc, :) = [];
            idx = 1;
        end
        point = point + 1;
    end
    %add last sequence to the orderedSequence array
%     fprintf('Final sequence was %d long\n and was inserted at %d \n',idx,seqStrtIdx)
    ordered_centerline = [ordered_centerline(1:seqStrtIdx-1,:); 
                          sequence(1:offShootIdx, :);
                          ordered_centerline(seqStrtIdx:filledTo,:)];   
end