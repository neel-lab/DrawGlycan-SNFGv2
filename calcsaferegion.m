function outline = calcsaferegion(mspos,safezone)
% CALCSAFEREGION: calculate the outline of the area occupied by the
% specified monosacs.
%
% Syntax:
% outline = calcsaferegion(mspos,safezone)
%
% Input:
% mspos: n x 2 numerical array, the monosaccharide positions.
% safezone: the serial number of the monosac. that form the area.
%
% Output:
% outline: n x 3 numerical array, 1st column contains y values,
% 2nd column contains the x values of the left boundary, 3rd column
% contains the x values of the right boundary. Y values increases from up
% to down.
%
% Note:
% 'outline' contains the position of the points that formed the
% outline of the area. If the area contains only 1 point, the perimeter
% will be this point.
% Due to program limitations, concaves in the outline can only be partially
% described: program is not able to describe concaves at the up/down side
% of the area, but those in left/right side can be shown. Voids inside the
% area is not available either.
%
% Example:
% N/A. Use breakpoints in parent program to test.
%
% Children function:
% N/A
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


if ~isempty(safezone)
    safezone = unique(safezone);
    allsafemspos = mspos(safezone,:);
    [uniy,~,uniyind] = unique(allsafemspos(:,2));  % pick all unique y values, this will decide the size of output
    outline = zeros(max(uniyind),3);
    outline(:,1) = uniy;
    for i = 1:max(uniyind)
        thisrow = allsafemspos(uniyind == i,1);
        outline(i,2) = min(thisrow);  % min. x value of the monosac that shares the y value is the left boundary
        outline(i,3) = max(thisrow);  % max. is the right boundary.
    end
else
    outline = [];
end
end