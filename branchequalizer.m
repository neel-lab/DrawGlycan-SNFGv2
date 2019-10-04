function balmspos = branchequalizer(mspos,bondmap)
% BRANCHEQUALIZER: evenly distribute branches in a tree structure.
%
% Syntax:
% balmspos = branchequalizer(mspos,bondmap)
%
% Input:
% mspos: n x 2 numerical array, the monosaccharide position (not
% equalized).
% bondmap: n x n array of 0's and 1's, the linkage map of the glycan.
%
% Output:
% balmspos: the equalized monosaccharide positions.
%
% Note:
% Do not leave [0,0]'s in array except for the first monosac. as [0,0] is
% used to mark the origin. Bone structure is an exception.
%
% Example:
% N/A
%
% Children function:
% GLYTREETRACKER, CALCSAFEREGION, ANTICOLLISION
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


balmspos = mspos;
%% squeeze the tree
[~,~,unixind] = unique(mspos(:,1));
for i = 1:max(unixind)  % maximum compression
    [~,yseq] = sort(balmspos(unixind == i,2));
    newyval = 0:sum(unixind == i) - 1;
    balmspos(unixind == i,2) = newyval(yseq);
end
counter = 2;
if length(balmspos) > 2
    for i = max(balmspos(:,1)):-1:2
        thiscolumn = find(balmspos(:,1) == i);
        whosechild = arrayfun(@(x) find(bondmap(:,x)),thiscolumn);
        [uniqueparent,~,uniparind] = unique(whosechild);
        for j = max(uniparind):-1:1
            %% alt. balance method 1
            thesebranches = thiscolumn(uniparind == j);
            if mod(length(thesebranches),2)
                balpos = median(balmspos(thesebranches,2));
            else
                balpos = mean(balmspos(thesebranches,2));
            end
            
            %% alt. bal method 2
            %             balpos = mean(balmspos(thiscolumn(uniparind == j),2));
            %% end alt. bal method
            
            yfix = balpos - balmspos(uniqueparent(j),2);
            balmspos(uniqueparent(j),2) = balmspos(uniqueparent(j),2) + yfix;
        end
        %% compress this column, if applicable
        adjustedcolumn = find(balmspos(:,1) == i - 1);
        [~,ind] = sort(mspos(adjustedcolumn,2));
        adjustedcolumn = adjustedcolumn(ind);
        localsafeind = zeros(size(mspos,1),1);
        for j = 1:length(adjustedcolumn)
            currentarm = glytreetracker(bondmap,adjustedcolumn(j),[],'down');
            if j > 1
                % raw adjustment, bring those too high down, too low up, so
                % every monosac in this column have +1 spacing, only after
                % this should the anticollision be carried out, otherwise
                % there could be no collision detected.
                yfix = balmspos(adjustedcolumn(j),2) - balmspos(adjustedcolumn(j - 1),2) - 1;
                balmspos(currentarm,2) = balmspos(currentarm,2) - yfix;
            end
            currentzone = calcsaferegion(balmspos,currentarm);
            safezone = calcsaferegion(balmspos,find(localsafeind));
            yfix = anticollision(safezone,currentzone,'up');
            balmspos(currentarm,2) = balmspos(currentarm,2) + yfix;
            localsafeind(currentarm) = 1;
        end
        counter = counter + 1;
    end
    %% handle stub
    %% alt. bal method 1
    thesebranches = find(bondmap(1,:) == 1);
    if mod(length(thesebranches),2)
        balmspos(1,2) = median(balmspos(thesebranches,2));
    else
        balmspos(1,2) = mean(balmspos(thesebranches,2));
    end
    
    %% alt bal method 2
    %     balmspos(1,2) = mean(balmspos(bondmap(1,:) == 1,2));
    %% end alt. bal method
end
origin = balmspos(1,:);
balmspos = balmspos - repmat(origin,size(balmspos,1),1);  % normalize all coordinates
end