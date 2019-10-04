function output = drawglycan(allsgp,varargin)
% DRAWGLYCAN: draw glycan or glycopeptide. This is the main program
%
% Syntax:
% output = drawglycan(allsgp,options)
% output = drawglycan(allsgp,optionname1,optionvalue1,...)
%
% Input:
% allsgp: string, glycan and glycopeptide string in IUPAC format.
% options: structure, user defined options listed below in Table format.
% optionname and optionvalue: user defined option name and value pairs.
%
% Output:
% output: structure. As of now, it has a field "tcp", an 2 x 2 numerical
% array. First row is the coordinate of figure's lower left corner, second
% row is its upper right corner.
%
%
% Available options (case sensitive):
% ________________________________________________________________________________________
% |  option name              purpose                   option value                     |
% |--------------------------------------------------------------------------------------|
% |  orientation           Orientation of the        String:'up','down','left','right'.  |
% |                        glycan structure          Default is 'right'.                 |
% |  monosacsize           Monosac. size             Number in range: [0,1].             |
%                                                    Default=0.5.                        |
% |  msperiwidth           Monosac. border           Any integer. Default=1.             |
% |                        thickness                                                     |
% |  bondwidth             Molecular bond            Any integer. Default=2.             |
% |                        thickness                                                     |
% |  perpendicularmonosac  List of perpendicular     1 x n cell array of strings.        |
% |                        monosacs.                 Default= {'Fuc','Xyl'}.             |
% |  showlink              Diplay linkage info.      String, 'yes' or 'no'.              |
% |                                                  Default = 'yes'.                    |
% |  fontsize              Linkage data font size.   Any integer. Default = 12.          |
% |  linkinfodist          Linkage text position.    Any number. Default = 0.7.          |
% |  linkinfotheta         Text orientation with     Any number in range = [-90,90].     |
% |                        respect to linkage.       Default = 30 (degrees).             |
% |  bondbreaksiglength    Length of glycosidic      Any number. Default = 0.5.          |
% |                        bond fragmentation line.                                      |
% |  structspacing         Spacing between glycans   Any number. Default = 1.            |
% |                        listed in cell array.                                         |
% |  aaspacing             Spacing between amino     Any number. Default = 0.75.         |
% |                        acids.                                                        |
% |  fileout               File saving specs.        String, output figure file          |
% |                                                  name. Default = empty.              |
% |  visible               Display figure.           String,'on' or 'off',Default='on'.  |
% ----------------------------------------------------------------------------------------
% Example:
% drawglycan('A[GalNAc(a1-3)](-N "B2" -C  "X4")AB[Xyl(b1-3)GalNAc(a1-3)]B(-C "X2.5")C[Fuc(a1-3)]CDD')
%
% Children function:
% USG, CALCGLYPOS, ESTIFIGSIZE, PAINTGLYCAN, GETFRAGINFO, REARRANGEGLYPEP,
% PLOTPEPBREAK
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


% v2.0 update note: additional local options added, see "DrawGlycanPara"
% for details.
% v1.1 update note: draw one structure at a time, the function of drawing
% multiple structure from a cell input has been cancelled.

%% INITIATING - BUILD DEFAULT OPTIONS
defoptname = {'orientation','monosacsize','msperiwidth',...
    'bondwidth','perpendicularmonosac',...
    'linkinfodist','linkinfotheta','bondbreaksiglength',...
    'showlink','fontsize','workingmode',...
    'structspacing','aaspacing','fileout',...
    'visible','inputformat','displaybrackettext',...
    'specialoptions','linkfontsize','sortbranches',...
    'figurehandle','resize','pointsperunit'};  % option names

defoptvalue = {'left',.5,1,...
    2,{'Fuc','Xyl','Deoxyhexose','Pentose'},...
    .7,-30,.5,...
    'yes',8,'',...
    1,.75,'',...
    'on','iupac','no',...
    {},16,true,...
    [],true,50};  % option default values

defopt = usg('gen',defoptname,defoptvalue);  % defopt is the option structure, all values are now default
if ~isempty(varargin)
    if isstruct(varargin{1})
        % update option values with user input (input is a structure too)
        options = usg('mod',defopt,varargin{1});
    elseif ischar(varargin{1})
        % update option values (input is "option name" - "option value" pairs)
        options = usg('mod',defopt,varargin(1:2:end),varargin(2:2:end));
    elseif iscell(varargin{1})
        % update option values with user input (input is a cell array: {option1,value1,option2,value2,...})
        inputoption = varargin{1};
        options = usg('mod',defopt,inputoption(1:2:end),inputoption(2:2:end));
    end
else
    options = defopt;
end

if ismember(lower(options.orientation),{'right','down'})
    % program is designed around "up" and "left", here the angles need to
    % be processed to avoid mirroring
    options.linkinfotheta = -options.linkinfotheta;
end
options.isaxes = true;  % assuming we have the whole window for us to draw on
% otherwise it's a small area in the window

[glyseq,pepseq,glypos,PTMid] = distggp(allsgp,options.inputformat);  % split glycan and peptide from input sequence
% glyseq is cell, pepseq is char
if isempty(options.workingmode)
    if isempty(pepseq)
        options.workingmode = 'G';
    else
        if ~isempty(glyseq)
            options.workingmode = 'GP';
        else
            options.workingmode = 'P';
        end
    end
end
%% handle different input format, produce a "glysgp" and "specialoptions"
glyoptions = [DrawGlycanPara.intglymodinfo,DrawGlycanPara.intglybondinfo,...
    DrawGlycanPara.extglymodinfo,DrawGlycanPara.glyidentityinfo,'CURLYBRACKET'];

if strcmpi(options.inputformat,'IUPAC')
    switch options.workingmode
        case 'G'
            [cleanseq,specialoptions] = getglydressing(strtrim(glyseq));
            glysgp = IUPAC2Sgp(cleanseq,2);
        case 'GP'
            glyseq = cellfun(@strtrim,glyseq,'UniformOutput',false);
            [cleanseq,specialoptions] = cellfun(@getglydressing,glyseq,'uniformoutput',false);
            glysgp = cellfun(@IUPAC2Sgp,cleanseq,num2cell(repmat(2,numel(glyseq),1)),'uniformoutput',false);
        case 'P'
            specialoptions = options.specialoptions;
    end
    if ~isempty(specialoptions)  % Program is designed based on SGP format. The monosac. 
        % are arranged differently. Therefore, option's position needs to be adjusted
        if iscell(specialoptions)
            % This means there are multiple glycans to be drawn
            for i = 1:length(specialoptions)
                if ~isempty(specialoptions{i})
                    spoptfld = fieldnames(specialoptions{i});  % special option fields
                    numMono = length(strfind(glysgp{i},'{'));
                    for j = 1:length(spoptfld)
                        if ismember(upper(spoptfld{j}),glyoptions)
                            tempspopt = specialoptions{i}.(spoptfld{j});
                            for k = 1:size(tempspopt,1)
                                tempspopt{k,2} = bsxfun(@minus,numMono + 1,tempspopt{k,2});  % convert option associated monosac location from IUPAC to SGP
                            end
                            specialoptions{i}.(spoptfld{j}) = tempspopt;
                        end
                    end
                end
            end
        elseif isstruct(specialoptions)
            % Only 1 glycan
            spoptfld = fieldnames(specialoptions);  % special option fields
            numMono = length(strfind(glysgp,'{'));
            for j = 1:length(spoptfld)
                if ismember(upper(spoptfld),glyoptions)
                    tempspopt = specialoptions.(spoptfld{j});
                    for k = 1:size(tempspopt,1)
                        tempspopt{k,2} = bsxfun(@minus,numMono + 1,tempspopt{k,2});  % convert option associated monosac location from IUPAC to SGP
                    end
                    specialoptions.(spoptfld{j}) = tempspopt;
                end
            end
        end
    end
elseif strcmpi(options.inputformat,'GLYCAM')
    % NOT FINISHED YET - STILL IN TEST PHASE
    if ~any(strfind(allsgp,'('))
        allsgp = glycam2iupac(allsgp);  % not tested
    end
    options.inputformat = 'IUPAC';
    [cleanseq,specialoptions] = getglydressing(strtrim(allsgp));
    glysgp = IUPAC2Sgp(cleanseq,2);
elseif strcmpi(options.inputformat,'SGP1')
    switch options.workingmode
        case 'G'
            glysgp = sgp1to2(allsgp,PTMid);
            specialoptions = options.specialoptions;
        case 'GP'
            glyseq = cellfun(@strtrim,glyseq,'UniformOutput',false);
            glysgp = cellfun(@sgp1to2,glyseq,num2cell(PTMid),'UniformOutput',false);
            specialoptions = options.specialoptions;
            if isempty(specialoptions)
                specialoptions = cell(size(glysgp));
            end
        case 'P'
            specialoptions = options.specialoptions;
    end
elseif strcmpi(options.inputformat,'SGP2')
    switch options.workingmode
        case 'G'
            glysgp = getglydressing(strtrim(glyseq));
            specialoptions = options.specialoptions;
        case 'GP'
            glysgp = cellfun(@strtrim,glyseq,'UniformOutput',false);
            specialoptions = options.specialoptions;
            if isempty(specialoptions)
                specialoptions = cell(size(glysgp));
            end
    end
elseif strcmpi(options.inputformat,'LINUCS')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'BCSDB')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'LINEARCODE')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'GLYCOCT')
    % TO BE COMPLETED
    % sequence itself is a connection table, reformat and put it in specialoptions.
elseif strcmpi(options.inputformat,'WURCS')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'KCF')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'CABOSML')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'GLYDE')
    
else
    return
end
switch upper(options.workingmode)
    case 'G'
        [plotinfoout,specialoptions] = calcglypos(glysgp,options,specialoptions);  % only draw 1 glycan structure, for multiple structures, call func. multiple times
        % specialoptions are merged into plotinfoout
        if ~isempty(specialoptions)
            plotinfoout = calcattachmentpos(plotinfoout,options);
        end
        if isempty(options.figurehandle)
            parentfigure = figure;
            options.figurehandle = axes(parentfigure,'Units','normalized','Position',[0 0 1 1]);
            hold on;
            set(options.figurehandle,'visible',options.visible)
            options.isaxes = false;
        end
        fig = options.figurehandle;
        [plotinfoout,glytcp] = paintglycan(plotinfoout,options,[]);
        atmtcp = paintattachment(plotinfoout,options,glytcp);
        output.tcp = [min(glytcp(1,1),atmtcp(1,1)),min(glytcp(1,2),atmtcp(1,2));...
            max(glytcp(2,1),atmtcp(2,1)),max(glytcp(2,2),atmtcp(2,2))];
        estimatefigsize(plotinfoout,'',output.tcp,options);
    case 'GP'  % default format: IUPAC
        options.orientation = 'up';
        plotinfoout = cell(size(glysgp));
        for i = 1:length(glysgp)
            [plotinfoout{i},specialoptions{i}] = calcglypos(glysgp{i},options,specialoptions{i});  % only draw 1 glycan structure, for multiple structures, call func. multiple times
            % specialoptions are merged into plotinfoout
            if ~isempty(specialoptions{i})
                plotinfoout{i} = calcattachmentpos(plotinfoout{i},options);
            end
        end
        [pepoptions,cleanpepseq] = getseqoptions(pepseq,'p');
        if ~isempty(specialoptions{1}) && ismember('pepC',fieldnames(specialoptions{1}))
            if isempty(pepoptions) || ~isfield(pepoptions,'C')
                pepoptions.C = specialoptions{1}.pepC;
            else
                pepoptions.C = [pepoptions.C;specialoptions{1}.pepC];
            end
        end
        if ~isempty(specialoptions{1}) && ismember('pepN',fieldnames(specialoptions{1}))
            if isempty(pepoptions) || ~isfield(pepoptions,'N')
                pepoptions.N = specialoptions{1}.pepN;
            else
                pepoptions.N = [pepoptions.N;specialoptions{1}.pepN];
            end
        end
        
        [plotinfoout,peppos] = rearrangeglypep(cleanpepseq,glypos,plotinfoout,options);
        if isempty(options.figurehandle)
            parentfigure = figure;
            options.figurehandle = axes(parentfigure,'Units','normalized','Position',[0 0 1 1]);
            hold on;
            set(options.figurehandle,'visible',options.visible,'Units','normalized')
            options.isaxes = false;
        end
        fig = options.figurehandle;
        for i = 1:length(cleanpepseq)
            text(fig,peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        tcp = zeros(2,2);
        tcp(1,1) = peppos(1) - options.aaspacing;
        tcp(1,2) = -2;
        tcp(2,1) = peppos(end) + options.aaspacing;
        glytcp = tcp;
        atmtcp = tcp;
        for i = 1:length(plotinfoout)
            if ~isempty(plotinfoout{i})
                [plotinfoout{i},glytcp] = paintglycan(plotinfoout{i},options,glytcp);
                atmtcp = paintattachment(plotinfoout{i},options,atmtcp);
            end
        end
    case 'P'
        [pepoptions,cleanpepseq] = getseqoptions(pepseq,'p');
        peppos = options.aaspacing * [0:length(cleanpepseq)-1];
        tcp = zeros(2);
        tcp(1,1) = peppos(1) - options.aaspacing;
        tcp(1,2) = -2;
        tcp(2,1) = peppos(end) + options.aaspacing;
        glytcp = zeros(2);
        atmtcp = zeros(2);
        if isempty(options.figurehandle)
            parentfigure = figure;
            options.figurehandle = axes(parentfigure,'Units','normalized','Position',[0 0 1 1]);
            hold on;
            set(options.figurehandle,'visible',options.visible,'Units','normalized')
            options.isaxes = false;
        end
        fig = options.figurehandle;
        if ~isempty(specialoptions) && ismember('pepC',fieldnames(specialoptions{1}))
            pepoptions.C = specialoptions{1}.pepC;
        end
        if  ~isempty(specialoptions) && ismember('pepN',fieldnames(specialoptions{1}))
            pepoptions.N = specialoptions{1}.pepN;
        end
        for i = 1:length(cleanpepseq)
            text(fig,peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        glytcp = zeros(2);
        atmtcp = zeros(2);
        tcp = zeros(2);
end
if ~options.isaxes
    set(fig,'position',[0,0,1,1],'visible','off')
end
axis equal
if strcmpi(options.workingmode,'gp') || strcmpi(options.workingmode,'p')
    pbtcp = plotpepbreak(peppos,pepoptions,glypos,options);  % this is the last step because frag marker length is related to figure size
    if isfield(pepoptions,'NGMOD')
        ngmods = pepoptions.NGMOD;
        for i = 1:size(ngmods,1)
            text(fig,peppos(ngmods{i,2}),-2,ngmods{i,1},'FontName','Helvetica','VerticalAlignment','bottom',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
    end
    tcp = [min([glytcp(1,1),atmtcp(1,1),tcp(1,1),pbtcp(1,1)]),min([glytcp(1,2),atmtcp(1,2),tcp(1,2),pbtcp(1,2)]);...
        max([glytcp(2,1),atmtcp(2,1),tcp(2,1),pbtcp(2,1)]),max([glytcp(2,2),atmtcp(2,2),tcp(2,2),pbtcp(2,2)])];
    output.tcp = tcp;
    if strcmpi(options.workingmode,'gp')
        set(fig,'XLim',[tcp(1,1),tcp(2,1)],'YLim',[tcp(1,2),tcp(2,2)+options.monosacsize]);
    elseif strcmpi(options.workingmode,'p')
        set(fig,'XLim',[peppos(1) - options.aaspacing,peppos(i) + options.aaspacing],'YLim',[-2.5,1]);
    end
end

if ~isempty(options.fileout)
    set(gcf, 'InvertHardcopy', 'off')
    saveas(gcf,options.fileout)
end
if strcmpi(options.visible,'off')
    p=get(options.figurehandle,'Parent');
    close(p);
end
end