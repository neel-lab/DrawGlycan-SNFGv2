function [plotinfoout,peppos] = rearrangeglypep(pepseq,glycosylpos,plotinfoout,options)
% REARRANGEGLYPEP: adjust the relative position of glycan(s) and peptide
% backbone in glycopeptide. Provides position of amino acids in the peptide
% backbone in the DrawGlycan canvas.
%
% Syntax:
% [plotinfoout,peppos] = rearrangeglypep(pepseq,glycosylpos,plotinfoout,options)
%
% Input:
% pepseq: string, amino acid sequence of peptide backbone of glycopeptide.
% glycosylpos: 1 x n numerical array, the number of the amino acid that has
% glycan attached to it.
% plotinfoout: m x n cell array, contains the plot information used in
% painting program, 1st column contains the positioning info for each
% monosaccharide in each glycan.
% options: struct, here 'structspacing' and 'fontsize' is used to estimate
% figure size and relative position of peptide and glycan(s).
%
% Output:
% plotinfoout: m x n cell array, the modified plot information after
% consideration of relative position of glycan(s) and peptide backbone.
% peppos: 1 x n numerical array, the position of each amino acid in the
% peptide backbone. The numbers are the x position of amino acid. They will
% be put at [x,-1] with "verticalalignment" set to "top".
%
% Note:
% Here AAs' letter size is set to be 75% the diameter of monosac. symbols.
% If not affected by glycan spacing requirements, the spacing will be 0.5
% unit. The actual spacing in figure between AA's will
% then be decided only by axis resolution, or to be exact, the figure size
% in pixels. It is not recommended to use fontsize too big as overlapping text
% or figures may ensue.
%
% Example:
% N/A. Use breakpoints in main program DRAWGLYPEP to test this code.
%
% Children function:
% N/A

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


% change "mspos" in "plotinfoout" to compensate for AA string
% Here we assume each letter in peptide backbone has 75% the size as
% monosac
AAspacing = options.aaspacing;

peppos = zeros(length(pepseq),1);
%% estimate spacing between glycans
glyposspacing = diff(glycosylpos);

if isempty(glyposspacing)  % 1 glycan attached
    tempstr = pepseq(1:glycosylpos);
    previouspos = AAspacing;
    for i = length(tempstr):-1:1
        peppos(i) = previouspos - AAspacing;
        previouspos = peppos(i);
    end
    tempstr = pepseq(glycosylpos+1:end);
    previouspos = 0;
    for i = 1:length(tempstr)
        peppos(glycosylpos + i) = previouspos + AAspacing;
        previouspos = peppos(glycosylpos + i);
    end
else  % >= 2 glycans attached
    rootpos = zeros(size(plotinfoout,1),1);
    orimspos = plotinfoout{1}.mspos;
    orimsright = max(orimspos(:,1));
    fldnames = lower(fieldnames(plotinfoout{1}));
    if ismember('cbitem',fldnames)
        cbplotinfos = plotinfoout{1}.CBITEM.cbplotinfosto;
        cbright = max(plotinfoout{1}.CBITEM.cbplotinfosto(1).mspos(:,1));
        for j = 2:size(cbplotinfos,1)
            cbitempos = plotinfoout{1}.CBITEM.cbplotinfosto(j).mspos;
            cbitemxspan = max(cbitempos(:,1)) - min(cbitempos(:,1));
            cbright = cbright + cbitemxspan + options.structspacing;
        end
    else
        cbright = [];
    end
    orirightlimit = max([orimsright,cbright]);
    lastLR = [0,orirightlimit];
    allLR = [lastLR;zeros(length(plotinfoout)-1,2)];
    for i = 2:size(plotinfoout,1)
        thismspos = plotinfoout{i}.mspos;
        thismsLR = [min(thismspos(:,1)),max(thismspos(:,1))];
        fldnames = lower(fieldnames(plotinfoout{i}));
        if ismember('cbitem',fldnames)
            cbplotinfos = plotinfoout{i}.CBITEM.cbplotinfosto;
            cbLR = [min(plotinfoout{i}.CBITEM.cbplotinfosto(1).mspos(:,1)),max(plotinfoout{i}.CBITEM.cbplotinfosto(1).mspos(:,1))];
            for j = 2:size(cbplotinfos,1)
                cbitempos = plotinfoout{i}.CBITEM.cbplotinfosto(j).mspos;
                cbitemxspan = max(cbitempos(:,1)) - min(cbitempos(:,1));
                cbLR(2) = cbLR(2) + cbitemxspan + options.structspacing;
            end
        else
            cbLR = [];
        end
        tempallLR = [thismsLR;cbLR];
        thisLR = [min(tempallLR(:,1));max(tempallLR(:,2))];
        allLR(i,:) = thisLR;
        rootpos(i) = rootpos(i-1) + options.structspacing + lastLR(2) + abs(thisLR(1));
        lastLR = thisLR;
        % calculate the min. distance between glycans: 1st glycan as origin
        % point, add the width of the right half, min. spacing between
        % glycans, the left half width of the next glycan to get it.
    end
    glycosylposdistpxl = glyposspacing * AAspacing;
    rootspacing = diff(rootpos);
    mode = rootspacing >= glycosylposdistpxl;
    %% start deciding positioning of amino acids
    %% use the 1st glycosyl pos as anchor point [0,0]
    %% before it:
    tempstr = pepseq(1:glycosylpos(1));
    previouspos = AAspacing;
    for i = length(tempstr):-1:1
        peppos(i) = previouspos - AAspacing;
        previouspos = peppos(i);
    end
    %% for the rest:
    previouspos = 0;
    for i = 1:length(mode)
        if mode(i)  % true means need to stretch peptidebackbone
            tempstr = pepseq(glycosylpos(i)+1:glycosylpos(i+1));
            startend = [previouspos,previouspos + rootspacing(i)];
            localAAspacing = diff(startend)/length(tempstr);
            previouspos = startend(1);
            for j = 1:length(tempstr)
                previouspos = previouspos + localAAspacing;
                peppos(glycosylpos(i) + j) = previouspos;
            end
        else  % readjust glycan position according to AA position
            tempstr = pepseq(glycosylpos(i) + 1:glycosylpos(i+1));
            startend = [previouspos,previouspos+(glycosylpos(i+1)-glycosylpos(i))*AAspacing];
            rootpos(i+1) = rootpos(i) + diff(startend);
            previouspos = peppos(glycosylpos(i));
            for j = 1:length(tempstr)
                previouspos = previouspos + AAspacing;
                peppos(glycosylpos(i) + j) = previouspos;
            end
        end
    end
    rootpos = peppos(glycosylpos);
    for i = 2:length(plotinfoout)
        fldnames = lower(fieldnames(plotinfoout{i}));
        plotinfoout{i}.mspos = plotinfoout{i}.mspos + repmat([rootpos(i),0],size(plotinfoout{i}.mspos,1),1);
        if ismember('cbitem',fldnames)
            plotinfoout{i}.CBITEM.bracketwaypoints = plotinfoout{i}.CBITEM.bracketwaypoints + repmat([rootpos(i),0],size(plotinfoout{i}.CBITEM.bracketwaypoints,1),1);
            plotinfoout{i}.CBITEM.cbplotinfosto(1).mspos = plotinfoout{i}.CBITEM.cbplotinfosto(1).mspos + repmat([rootpos(i),0],size(plotinfoout{i}.CBITEM.cbplotinfosto(1).mspos,1),1);
        end
    end
    %% handle tail, if any
    if length(pepseq) > glycosylpos(end)
        tempstr = pepseq(glycosylpos(end)+1:end);
        previouspos = peppos(glycosylpos(end));
        for i = 1:length(tempstr)
            previouspos = previouspos + AAspacing;
            peppos(glycosylpos(end) + i) = previouspos;
        end
    end
end
% ADDUCT - COMBINE & REPOSITION
alladducttext = '';
adductpos = [];
bracketpos = [];
for i = 1:length(plotinfoout)
    if isfield(plotinfoout{i},'ADDUCT')
        for j = 1:size(plotinfoout{i}.ADDUCT,1)
            alladducttext = [alladducttext,plotinfoout{i}.ADDUCT{j,1},' '];
        end
        adductpos = [adductpos;i];
        bracketpos = [bracketpos;plotinfoout{i}.bracketspos];
    end
end
if ~isempty(alladducttext) && strcmpi(alladducttext(end),' ')
    alladducttext = alladducttext(1:end-1);
    brackettop = max(bracketpos(:,2));
    newbracketpos = [peppos(1) - options.aaspacing + .1,brackettop;...
        peppos(1) - options.aaspacing,brackettop;...
        peppos(1) - options.aaspacing,-1;...
        peppos(1) - options.aaspacing + .1,-1;...
        peppos(end) + options.aaspacing - .1,brackettop;...
        peppos(end) + options.aaspacing,brackettop;...
        peppos(end) + options.aaspacing,-1;...
        peppos(end) + options.aaspacing - .1,-1];
    newadducttextpos = [peppos(end) + options.aaspacing + .1,brackettop/2-.5];
    for i = 1:length(adductpos)
        if i == 1
            plotinfoout{adductpos(i)}.bracketspos = newbracketpos;
            plotinfoout{adductpos(i)}.adducttextpos = newadducttextpos;
            plotinfoout{adductpos(i)}.ADDUCT = {alladducttext,1};
        else
            plotinfoout{adductpos(i)} = rmfield(plotinfoout{adductpos(i)},{'ADDUCT','bracketspos','adducttextpos'});
        end
    end
end