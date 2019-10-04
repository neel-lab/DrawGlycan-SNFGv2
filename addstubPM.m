function [mspos,directionseq] = addstubPM(mspos,bondmap,stub,bonebondmap,PMroot,pmposstorage,pmindstorage,safezone,directionseq)
% ADDSTUBPM: add perpendicular subtree(s) to the stub, i.e. the region from
% the first monosac. to the first branching point in the sequence.
%
% Syntax:
% [mspos,directionseq] = addstubPM(mspos,bondmap,stub,bonebondmap,PMroot,pmsposstorage,pmsindstorage,safezone,directionseq)
%
% Input:
% mspos: n x 2 numerical array, the position of monsac. in glycan.
% bondmap: n x n numerical array, the linkage map of monosac. in glycan.
% stub
% bonebondmap: n x n numerical array, the linkage map of bone glycan
% structure. Has the same size as 'bondmap' but contains limited info.
% PMroot: m x 1 numerical array, the serial number of monosac. that has PM
% attached.
% pmposstorage: m x 2 cell array of k x 2 numerical array, the position
% info of PM subtrees.
% pmindstorage: m x 2 cell array of k x 1 numerical array, the serial
% number of monosac in PM subtrees.
% safezone: n x 1 numerical array, indicator of whether the position of the
% corresponding monosac. has been established and fixed.
% directionseq: 1 x n numerical array, the orientation of monosac. symbols.
%
% Output:
% mspos: modified position info of monosac. in glycan after combination of
% PM subtree and bone structure
% directionseq: modified orientation info.
%
% Note:
% N/A
%
% Example:
% N/A. Set breakpoints in main program.
%
% Children function:
% GLYTREETRACKER, CALCSAFEREGION, ANTICOLLISION
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


if stub ~= 0
    allparent = glytreetracker(bondmap,stub,[],'up');
    allparent = [allparent;stub];
    allparent = intersect(allparent,PMroot);
    for i = length(allparent):-1:1
        thisms = allparent(i);
        haschild = bonebondmap(thisms,:) ~= 0;  % does this monosac has other ms attached to it
        if any(haschild)
            mspos(thisms,2) = mean(mspos(haschild,2));
        else
            mspos(thisms,2) = 0;  % if this is the end of linear struct, give it a 0 as y-value
        end
        safezone(thisms) = 1;
    end
    if any(ismember(PMroot,allparent))
        stubroot = find(ismember(PMroot,allparent));
        for i = length(stubroot):-1:1
            instub = PMroot(stubroot(i));
            tempPMSstruct = pmposstorage(stubroot(i),:);
            if ~isempty(tempPMSstruct{1})  % if there is a upper PMS subtree
                upperpmpos = tempPMSstruct{1}+ repmat(mspos(instub,:),size(tempPMSstruct{1},1),1);
                uppersafezone = calcsaferegion(upperpmpos,1:size(upperpmpos,1));
            else
                upperpmpos = [];
                uppersafezone = [0 0 0];
            end
            if ~isempty(upperpmpos)
                exsafezone = calcsaferegion(mspos,find(safezone));
                tempexsafezone = [-exsafezone(:,1),exsafezone(:,2:3)];
                tempuppersafezone = [-uppersafezone(:,1),uppersafezone(:,2:3)];
                fix = anticollision(tempexsafezone,tempuppersafezone,'up');
                upperpmpos(:,2) = upperpmpos(:,2) - fix;  % this "flipping" technique will make downwards adjustments.
                mspos(pmindstorage{stubroot(i),1},:) = upperpmpos;
            end
            safezone(pmindstorage{stubroot(i),1}) = 1;
            if ~isempty(tempPMSstruct{2})  % if there is a lower PMS subtree
                lowerpmpos = tempPMSstruct{2}+ repmat(mspos(instub,:),size(tempPMSstruct{2},1),1);
                lowersafezone = calcsaferegion(lowerpmpos,1:size(lowerpmpos,1));
            else
                lowerpmpos = [];
                lowersafezone = [0 0 0];
            end
            if ~isempty(lowerpmpos)
                exsafezone = calcsaferegion(mspos,find(safezone));
                fix = anticollision(exsafezone,lowersafezone,'up');
                lowerpmpos(:,2) = lowerpmpos(:,2) + fix;
                mspos(pmindstorage{stubroot(i),2},:) = lowerpmpos;
            end
            safezone(pmindstorage{stubroot(i),2}) = 1;
        end
    end
end