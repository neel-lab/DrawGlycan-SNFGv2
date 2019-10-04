function monosacpos = drawrawtree(bondmap,distance,skiplist)
% DRAWRAWTREE: generate the initial monosaccharide positions of a glycan in
% the form of a tree structure
%
% Syntax:
% monosacpos = drawrawtree(bondmap,distance,skiplist)
%
% Input:
% bondmap: n x n matrix, the linkage map of glycan, n equals to the number
% of monosaccharides in the glycan. Element (M,N) marks the existence of
% glycosidic bond between Mth and Nth monosaccharide in the SGP2.0 sequence
% of the glycan, 1 means there exists a bond, 0 means no.
% distance: n x 1 numerical array. The distance of the monosaccharide from
% the first glycan in the sequence.
% skiplist: m x 1 numerical array, the serial number of the
% monosaccharide(s) not to be included in the position calculation process.
%
% Output:
% monosacpos: n x 2 numerical array. Calculated monosaccharide position.
%
% Note:
% DRAWRAWTREE creates the initial guess of monosaccharide position, needs to
% be further processed to create the final structure.
% The horizontal position of the monosaccharide is purely decided by its
% 'distance' value, the vertical position is decided by which branch it
% belongs to.
%
% Example:
% N/A. Set breakpoints in parent program.
%
% Children function:
% N/A
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


monosacpos = zeros(length(distance),2);
vert = 0;
for i = 1:length(distance)
    if ~ismember(i,skiplist)
        isend = ~any(bondmap(i,:));  % zeros in bondmap means
        % there is no other monosac connected to the current one,
        % which means a branch is finished, need to start a new one.
        monosacpos(i,:) = [distance(i) - 1,vert];
        if isend
            vert = vert + 1;
        end
    end
end
end