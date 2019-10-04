function structout = usg(mode,varargin)
% USG: Universal option Structure Generator. This program collates the
% options used by DrawGlycan (or other software).
%
% Syntax:
% structout = usg('gen',parametername,defaultvalue)
% structout = usg('gen',name1,value1,name2,value2,...)
% structout = usg('mod',originalstruct,modstruct)
% structout = usg('mod',originalstruct,name1,value1,...)
%
% Input:
% 'gen': generates parameter structure with default values
% 'mod': modifies existing structure use information from a user input
% structure
% parametername: string or 1-D cell array of strings containing parameter names
% defaultvalue: value or 1-D cell array of values with correspondence with
% paramtername
% name1,name2,...: string, struct element name
% value1,value2,...: corresponding values.
% originalstruct: structure, containing the original structure.
% modstruct: structure, if fieldnames in this structure also exists in
% originalstruct, the corresponding values will be replaced by the ones in
% modstruct.
%
% Output:
% structout: structure with all options and corresponding values
%
% Note:
% When used inside other functions, if 'varargin' is used as option
% modifier container, 'varargin' can be used the same way as 'modstruct' as
% mentioned in Syntax.
%
% Example:
% usg('gen','a',6,'b',9,'c','h')
%
% ans =
%
%     a: 6
%     b: 9
%     c: 'h'
%
% usg('gen',{'a','b'},{zeros(3),'mk'})
%
% ans =
%
%     a: [3x3 double]
%     b: 'mk'
%
% a=usg('gen','a',6,'b',9,'c','h');
% b=usg('gen',{'a','b'},{zeros(3),'mk'});
% usg('mod',a,b)
%
% ans =
%
%     a: [3x3 double]
%     b: 'mk'
%     c: 'h'
%
%
% usg('mod',b,a)
%
% ans =
%
%     a: 6
%     b: 9
%
% Children function:
% N/A
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


if strcmpi(mode,'gen')
    if nargin == 3              % input could be string or cell
        if ischar(varargin{1})  % if it is string, change it to a cell
            varargin{1} = varargin(1);
        end
        if ~iscell(varargin{2})
            varargin{2} = varargin(2);
        end
        parametername = varargin{1};    % creates a cell
        defaultvalue = varargin{2};     % creates a cell
    else                        % this is for parts of inputs
        if mod(nargin,2)  == 0
            error('Check input, number of input not match')
        else
            numofpairs = (nargin - 1)/2;
            parametername = cell(numofpairs,1);     % preallocate memory
            defaultvalue = cell(numofpairs,1);
            for i = 1:numofpairs
                parametername{i} = varargin{i*2-1}; % put values into cell
                defaultvalue{i} = varargin{i*2};
            end
        end
    end
    if numel(parametername) ~= numel(defaultvalue)
        error('Check input, input not in pairs')
    else
        for i = 1:numel(parametername)
            structout.(parametername{i}) = defaultvalue{i}; % convert cell to structure
        end
    end
elseif strcmpi(mode,'mod')
    if ~isa(varargin{1},'struct')
        error('Input 1 is not a structure')
    end
    if isempty(varargin{2})
        error('Cannot handle empty modification.');
    end
    if isstruct(varargin{2})
        oristruct = varargin{1};
        modstruct = varargin{2};
        orifieldnames = fieldnames(oristruct);
        modfieldnames = fieldnames(modstruct);
        fieldtobemod = orifieldnames(ismember(orifieldnames,modfieldnames));
        for i = 1:numel(fieldtobemod)
            oristruct.(fieldtobemod{i}) = modstruct.(fieldtobemod{i});
        end
    elseif ischar(varargin{2})
        if mod(nargin,2)  ~= 0
            error('Check input, input not in pairs')
        end
        oristruct = varargin{1};
        modfieldnames = varargin(2:2:end);
        modvalues = varargin(3:2:end);
        orifieldnames = fieldnames(oristruct);
        modind = ismember(modfieldnames,orifieldnames);
        modfieldnames = modfieldnames(modind);
        modvalues = modvalues(modind);
        for i = 1:sum(modind)
            oristruct.(modfieldnames{i}) = modvalues{i};
        end
    elseif iscell(varargin{2})      % when used inside other function and input is 'varargin'
        oristruct = varargin{1};
        modfieldnames = varargin{2};
        modvalues = varargin{3};
        orifieldnames = fieldnames(oristruct);
        modind = ismember(modfieldnames,orifieldnames);
        modfieldnames = modfieldnames(modind);
        modvalues = modvalues(modind);
        for i = 1:sum(modind)
            oristruct.(modfieldnames{i}) = modvalues{i};
        end
    end
    structout = oristruct;
end
end
