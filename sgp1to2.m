function glysgp = sgp1to2(allsgp,PTMid)
% SGP1TO2: convert SGP 1.0 format sequence to SGP 2.0 standard
%
% Syntax:
% glysgp = sgp1to2(allsgp,PTMid)
%
% Input:
% allsgp: string, glycan sequence in SGP 1.0 format
% PTMid: number, 1 means this sequence is a glycan, 0 means it is a
% non-glycan PTM
%
% Output:
% glysgp: glycan sequence in SGP 2.0 format. Monosaccharides are converted
% using their generic names, e.g. Hexose, HexNAc,...
%
% Note:
% If input is non-glycan PTM (PTMid == 0 & starts with "<" & ends with ">"),
% the name will be copied without conversion.
%
% Example:
% 1. Glycan
% glysgp = sgp1to2('{n{h{s}}}',1)
%
% glysgp =
%
%     '{(??-?)HexNAc{(??-?)Hexose{(??-?)Deoxynonulosonate}}}'
%
% 2. Non-glycan
% glysgp = sgp1to2('<oxidation>',0)
%
% glysgp =
%
%     '{(??-?)<oxidation>}'
%
% Children function:
% N/A

% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved


sgp1name = {'n','h','s',...
    'f','x','g'...
    ,'u','k'};
sgp2name = {'(??-?)HexNAc','(??-?)Hexose','(??-?)Deoxynonulosonate',...
    '(??-?)Deoxyhexose','(??-?)Pentose','(??-?)Nonulosonate',...
    '(??-?)Hexuronate','(??-?)Nonulosonate'};
glysgp = '';
if PTMid == 1
    for i = length(allsgp):-1:1
        if ~ismember(allsgp(i),sgp1name)
            glysgp = [allsgp(i),glysgp];
        else
            for j = 1:length(sgp1name)
                if strcmpi(sgp1name{j},allsgp(i))
                    glysgp = [sgp2name{j},glysgp];
                    break
                end
            end
        end
    end
elseif PTMid == 0
    if strcmpi(allsgp(1),'<') && strcmpi(allsgp(end),'>')
        glysgp = ['{(??-?)',allsgp,'}'];
    end
end