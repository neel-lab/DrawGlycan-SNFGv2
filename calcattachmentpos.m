function plotinfoout = calcattachmentpos(plotinfoout,options)
% CALCATTACHMENTPOS: calculate position of sub-structures outside the main
% glycan structure.
%
% Syntax:
% plotinfoout = calcattachmentpos(plotinfoout,options)
%
% Input:
% plotinfoout: structure, the result from CALCGLYPOS. The calculated
% positioning info of glycan main structure and other attachments.
% options: structure, contains global options for drawing glycans. See
% DRAWGLYCAN for details.
%
% Output:
% plotinfoout: structure, with substructure info added. If applicable, new
% fields will be added:
% "CBITEM": structure, info of "fuzzy" structures. It has 2 fields:
% "cbplotinfosto": cell, each contain a structure describing one of the
% "fuzzy" part, this structure is similar to "plotinfoout".
% "breacketwaypoints": waypoints for drawing the curly bracket with lines
% and curves.
% "bracketspos": coordinate of square brackets (for displaying adduct ions)
% "adducttextpos": position of adduct ion text.
%
% Note:
% N/A
%
% Example:
% N/A. Run examples in DRAWGLYCAN and set breakpoints.
%
% Children function:
% N/A

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


addonfix = zeros(1,4);
%% 1. anomer
if isfield(plotinfoout,'ANOMER')
    anomer = plotinfoout.ANOMER;
    switch lower(options.orientation)
        case 'up'
            plotinfoout.mspos = [[0,0];bsxfun(@plus,plotinfoout.mspos,[0,1])];
        case 'down'
            plotinfoout.mspos = [[0,0];bsxfun(@plus,plotinfoout.mspos,[0,-1])];
        case 'left'
            plotinfoout.mspos = [[0,0];bsxfun(@plus,plotinfoout.mspos,[-1,0])];
        case 'right'
            plotinfoout.mspos = [[0,0];bsxfun(@plus,plotinfoout.mspos,[1,0])];
    end
    revbondmap = plotinfoout.bondmap;
    revbondmap = [[1,zeros(1,size(revbondmap,2)-1)];revbondmap];
    revbondmap = [zeros(size(revbondmap,1),1),revbondmap];
    plotinfoout.bondmap = revbondmap;
    plotinfoout.alllinkout = ['??-?';plotinfoout.alllinkout];
    plotinfoout.allms = [anomer{1};plotinfoout.allms];
    plotinfoout.directionseq = [0;plotinfoout.directionseq];
    bonddescription = plotinfoout.bonddescription;
    if ~isempty(bonddescription)
        bonddescription(:,1) = num2cell(cell2mat(bonddescription(:,1)) + 1);
    end
    plotinfoout.bonddescription = bonddescription;
    specialoptionnames = [DrawGlycanPara.intglybondinfo,'CHAR','RE','RS'];
    fldnames = fieldnames(plotinfoout);
    fldnames = fldnames(ismember(fldnames,specialoptionnames));
    for i = 1:length(fldnames)
        tempopt = plotinfoout.(fldnames{i});
        tempopt(:,2) = num2cell(cell2mat(tempopt(:,2)) + 1);
        plotinfoout.(fldnames{i}) = tempopt;
    end
    if isfield(plotinfoout,'CURLYBRACKET')
        tempopt = plotinfoout.CURLYBRACKET;
        for i = 1:size(tempopt,1)
            tempopt{i,2} = tempopt{i,2} + 1;
        end
        plotinfoout.CURLYBRACKET = tempopt;
    end
end

%% 2. curly brackets
if isfield(plotinfoout,'CURLYBRACKET')  % currently only support 1 curly bracket
    tr = DrawGlycanPara.tipradius;  % tip radius
    curlybracketitem = plotinfoout.CURLYBRACKET;
    mainmspos = plotinfoout.mspos;
    curlybracketcover = mainmspos(curlybracketitem{1,2},:);
    cbseq = curlybracketitem{1,1};
    if strcmp(cbseq(end),',')
        cbseq(end) = '';
    end
    [optcontent,optcontentstart,optcontentend] = regexp(cbseq,DrawGlycanPara.regexp_optionvalue,'match','start','end');
    if ~isempty(optcontent)
        for i = 1:length(optcontentstart)
            cbseq(optcontentstart(i):optcontentend(i)) = blanks(optcontentend(i) - optcontentstart(i) + 1);
        end
    end
    [optcontent,optcontentstart,optcontentend] = regexp(cbseq,DrawGlycanPara.regexp_monosaclinkage,'match','start','end');
    if ~isempty(optcontent)
        for i = 1:length(optcontentstart)
            cbseq(optcontentstart(i):optcontentend(i)) = blanks(optcontentend(i) - optcontentstart(i) + 1);
        end
    end
    commaind = strfind(cbseq,',');
    commaind = [[1;commaind(:) + 1],[commaind(:) - 1;length(cbseq)]];
    allcbseq = cell(1,size(commaind,1));
    for i = 1:size(commaind,1)
        allcbseq{i} = curlybracketitem{1,1}(commaind(i,1):commaind(i,2));
    end
    allcbseq = strtrim(allcbseq);
    cbplotinfosto = {};
    for i = 1:length(allcbseq)
        [cleanseq,specialoptions] = getglydressing(allcbseq{i});
        cbsgp = IUPAC2Sgp(cleanseq,2);
        if ~isempty(specialoptions)  % get option position in sgp seq, because it needs conversion
            spoptfld = fieldnames(specialoptions);  % special option fields
            numMono = length(strfind(cbsgp,'{'));
            for j = 1:length(spoptfld)
                tempspopt = specialoptions.(spoptfld{j});
                for k = 1:size(tempspopt,1)
                    tempspopt{k,2} = bsxfun(@minus,numMono + 1,tempspopt{k,2});  % convert option associated monosac location from IUPAC to SGP
                end
                specialoptions.(spoptfld{j}) = tempspopt;
            end
        end
        cbplotinfosto = [cbplotinfosto;calcglypos(cbsgp,options,specialoptions)];
    end
    structwidth = 0;
    
    switch lower(options.orientation)
        case 'up'
            cbspan = [min(curlybracketcover(:,1)),max(curlybracketcover(:,1))];
            for i = 1:length(cbplotinfosto)
                structwidth = structwidth + max(cbplotinfosto{i}.mspos(:,1)) - ...
                    min(cbplotinfosto{i}.mspos(:,1)) + options.structspacing;
            end
            horizontalfix =  (max(curlybracketcover(:,1)) + min(curlybracketcover(:,1)))/2 - ...
                (structwidth-options.structspacing)/2;
            for i = 1:length(cbplotinfosto)
                cbplotinfosto{i}.mspos(:,1) = cbplotinfosto{i}.mspos(:,1) + horizontalfix - ...
                    min(cbplotinfosto{i}.mspos(:,1) - cbplotinfosto{i}.mspos(1,1));
                cbplotinfosto{i}.mspos(:,2) = cbplotinfosto{i}.mspos(:,2) + max(curlybracketcover(:,2)) + 1;
                addonfix(1) = max([addonfix(1),max(cbplotinfosto{i}.mspos(:,2)-cbplotinfosto{i}.mspos(1,2))]);
            end
            addonfix(1) = addonfix(1) + 1;  % curly bracket takes 1 unit
            addonfix(2) = (max(curlybracketcover(:,1)) + min(curlybracketcover(:,1)))/2 + ...
                (structwidth-1)/2 - (max(mainmspos(:,1) - mainmspos(1,1)));
            addonfix(2) = addonfix(2) * (addonfix(2) > 0);
            addonfix(4) = (max(curlybracketcover(:,1)) + min(curlybracketcover(:,1)))/2 - ...
                (structwidth-1)/2 - (min(mainmspos(:,1) - mainmspos(1,1)));
            addonfix(4) = addonfix(4) * (addonfix(4) < 0);
            midpoint = mean([max(curlybracketcover(:,1)),min(curlybracketcover(:,1))]);
            if max(curlybracketcover(:,1)) == min(curlybracketcover(:,1))
                cbspan = [-.25 .25]+curlybracketcover(1,1);
            end
            bracketwaypoints = [min(cbspan),max(curlybracketcover(:,2)) + .5 - tr;...
                min(cbspan),max(curlybracketcover(:,2))+.5;...
                midpoint-tr,max(curlybracketcover(:,2)) + .5;...
                midpoint-tr,max(curlybracketcover(:,2)) + .5 + tr;...
                midpoint+tr,max(curlybracketcover(:,2)) + .5 + tr;...
                midpoint+tr,max(curlybracketcover(:,2)) + .5;...
                max(cbspan),max(curlybracketcover(:,2)) + .5;
                max(cbspan),max(curlybracketcover(:,2)) + .5 - tr];
        case 'down'
            cbspan = [min(curlybracketcover(:,1)),max(curlybracketcover(:,1))];
            for i = 1:length(cbplotinfosto)
                structwidth = structwidth + max(cbplotinfosto{i}.mspos(:,1)) - ...
                    min(cbplotinfosto{i}.mspos(:,1)) + options.structspacing;
            end
            horizontalfix = (max(curlybracketcover(:,1)) + min(curlybracketcover(:,1)))/2 - ...
                (structwidth-options.structspacing)/2;
            
            for i = 1:length(cbplotinfosto)
                cbplotinfosto{i}.mspos(:,1) = cbplotinfosto{i}.mspos(:,1) + horizontalfix - ...
                    min(cbplotinfosto{i}.mspos(:,1) - cbplotinfosto{i}.mspos(1,1));
                cbplotinfosto{i}.mspos(:,2) = cbplotinfosto{i}.mspos(:,2) + min(curlybracketcover(:,2)) - 1;
                addonfix(3) = min([addonfix(3),min(cbplotinfosto{i}.mspos(:,2)-cbplotinfosto{i}.mspos(1,2))]);
            end
            addonfix(3) = addonfix(3) - 1;
            addonfix(2) = (max(curlybracketcover(:,1)) + min(curlybracketcover(:,1)))/2 + ...
                (structwidth-1)/2 - (max(mainmspos(:,1) - mainmspos(1,1)));
            addonfix(2) = addonfix(2) * (addonfix(2) > 0);
            addonfix(4) = (max(curlybracketcover(:,1)) + min(curlybracketcover(:,1)))/2 - ...
                (structwidth-1)/2 - (min(mainmspos(:,1) - mainmspos(1,1)));
            addonfix(4) = addonfix(4) * (addonfix(4) < 0);
            midpoint = mean([max(curlybracketcover(:,1)),min(curlybracketcover(:,1))]);
            if max(curlybracketcover(:,1)) == min(curlybracketcover(:,1))
                cbspan = [-.25 .25]+curlybracketcover(1,1);
            end
            bracketwaypoints = [min(cbspan),min(curlybracketcover(:,2)) - .5 + tr;...
                min(cbspan),min(curlybracketcover(:,2)) - .5;...
                midpoint-tr,min(curlybracketcover(:,2)) - .5;...
                midpoint-tr,min(curlybracketcover(:,2)) - .5 - tr;...
                midpoint+tr,min(curlybracketcover(:,2)) - .5 - tr;...
                midpoint+tr,min(curlybracketcover(:,2)) - .5;...
                max(cbspan),min(curlybracketcover(:,2)) - .5;
                max(cbspan),min(curlybracketcover(:,2)) - .5 + tr];
            
        case 'left'
            cbspan = [min(curlybracketcover(:,2)),max(curlybracketcover(:,2))];
            for i = 1:length(cbplotinfosto)
                structwidth = structwidth + max(cbplotinfosto{i}.mspos(:,2)) - ...
                    min(cbplotinfosto{i}.mspos(:,2)) + options.structspacing;
            end
            verticalfix = (max(curlybracketcover(:,2)) + min(curlybracketcover(:,2)))/2 - ...
                (structwidth-options.structspacing)/2;
            cbplotinfosto{1}.mspos(:,2) = cbplotinfosto{1}.mspos(:,2) + verticalfix - ...
                min(cbplotinfosto{1}.mspos(:,2) - cbplotinfosto{1}.mspos(1,2));
            for i = 1:length(cbplotinfosto)
                cbplotinfosto{i}.mspos(:,1) = cbplotinfosto{i}.mspos(:,1) + min(curlybracketcover(:,1)) - 1;
                addonfix(4) = min([addonfix(4),min(cbplotinfosto{i}.mspos(:,1)-cbplotinfosto{i}.mspos(1,1))]);
            end
            addonfix(4) = addonfix(4) - 1;
            addonfix(1) = (max(curlybracketcover(:,2)) + min(curlybracketcover(:,2)))/2 + ...
                (structwidth-1)/2 - (max(mainmspos(:,2) - mainmspos(1,2)));
            addonfix(1) = addonfix(1) * (addonfix(1) > 0);
            addonfix(3) = (max(curlybracketcover(:,2)) + min(curlybracketcover(:,2)))/2 - ...
                (structwidth-1)/2 - (min(mainmspos(:,2) - mainmspos(1,2)));
            addonfix(3) = addonfix(3) * (addonfix(3) < 0);
            midpoint = mean([max(curlybracketcover(:,2)),min(curlybracketcover(:,2))]);
            if max(curlybracketcover(:,2)) == min(curlybracketcover(:,2))
                cbspan = [-.25 .25]+curlybracketcover(1,2);
            end
            bracketwaypoints = [min(curlybracketcover(:,1)) - .5 + tr,min(cbspan);...
                min(curlybracketcover(:,1)) - .5,min(cbspan);...
                min(curlybracketcover(:,1)) - .5,midpoint - tr;...
                min(curlybracketcover(:,1)) - .5 - tr,midpoint - tr;...
                min(curlybracketcover(:,1)) - .5 - tr,midpoint + tr;...
                min(curlybracketcover(:,1)) - .5,midpoint + tr;...
                min(curlybracketcover(:,1)) - .5,max(cbspan);
                min(curlybracketcover(:,1)) - .5 + tr,max(cbspan)];
            
        case 'right'
            cbspan = [min(curlybracketcover(:,2)),max(curlybracketcover(:,2))];
            for i = 1:length(cbplotinfosto)
                structwidth = structwidth + max(cbplotinfosto{i}.mspos(:,2)) - ...
                    min(cbplotinfosto{i}.mspos(:,2)) + options.structspacing;
            end
            verticalfix = (max(curlybracketcover(:,2)) + min(curlybracketcover(:,2)))/2 - ...
                (structwidth-options.structspacing)/2;
            for i = 1:length(cbplotinfosto)
                cbplotinfosto{i}.mspos(:,2) = cbplotinfosto{i}.mspos(:,2) + verticalfix - ...
                    min(cbplotinfosto{i}.mspos(:,2) - cbplotinfosto{i}.mspos(1,2));
                cbplotinfosto{i}.mspos(:,1) = cbplotinfosto{i}.mspos(:,1) + max(curlybracketcover(:,1)) + 1;
                addonfix(2) = max([addonfix(2),max(cbplotinfosto{i}.mspos(:,1)-cbplotinfosto{i}.mspos(1,1))]);
            end
            addonfix(2) = addonfix(2) + 1;
            addonfix(1) = (max(curlybracketcover(:,2)) + min(curlybracketcover(:,2)))/2 + ...
                (structwidth-1)/2 - (max(mainmspos(:,2) - mainmspos(1,2)));
            addonfix(1) = addonfix(1) * (addonfix(1) > 0);
            addonfix(3) = (max(curlybracketcover(:,2)) + min(curlybracketcover(:,2)))/2 - ...
                (structwidth-1)/2 - (min(mainmspos(:,2) - mainmspos(1,2)));
            addonfix(3) = addonfix(3) * (addonfix(3) < 0);
            midpoint = mean([max(curlybracketcover(:,2)),min(curlybracketcover(:,2))]);
            if max(curlybracketcover(:,2)) == min(curlybracketcover(:,2))
                cbspan = [-.25 .25]+curlybracketcover(1,2);
            end
            bracketwaypoints = [max(curlybracketcover(:,1)) + .5 - tr,min(cbspan);...
                max(curlybracketcover(:,1)) + .5,min(cbspan);...
                max(curlybracketcover(:,1)) + .5,midpoint - tr;...
                max(curlybracketcover(:,1)) + .5 + tr,midpoint - tr;...
                max(curlybracketcover(:,1)) + .5 + tr,midpoint + tr;...
                max(curlybracketcover(:,1)) + .5,midpoint + tr;...
                max(curlybracketcover(:,1)) + .5,max(cbspan);
                max(curlybracketcover(:,1)) + .5 - tr,max(cbspan)];
    end
    fldnames = {};  % unifying structure
    for i = 1:length(cbplotinfosto)
        fldnames=[fldnames;fieldnames(cbplotinfosto{i})];
    end
    fldnames = unique(fldnames);
    combicbplotinfosto = [];
    for i = 1:length(cbplotinfosto)
        for j = 1:length(fldnames)
            if isfield(cbplotinfosto{i},fldnames{j})
                combicbplotinfosto(i).(fldnames{j}) = cbplotinfosto{i}.(fldnames{j});
            else
                combicbplotinfosto(i).(fldnames{j}) = [];
            end
        end
    end
    cbitem.cbplotinfosto = combicbplotinfosto;
    cbitem.bracketwaypoints = bracketwaypoints;
    plotinfoout.CBITEM = cbitem;
end

%% 3. adduct
if isfield(plotinfoout,'ADDUCT')
    mspos = plotinfoout.mspos;
    if min(mspos(:,2)) == max(mspos(:,2)) && addonfix(1) == 0 && addonfix(3) == 0  % for linear structure, if there are curly bracket items usually addonfix will handle it
        addonfix(1) = addonfix(1) + .5;
        addonfix(3) = addonfix(3) - .5;
    end
    bracketspos = zeros(8,2);
    bracketspos(1:4,1) = min(mspos(:,1)) + addonfix(4) - [.4;.5;.5;.4];
    bracketspos(5:8,1) = max(mspos(:,1)) + addonfix(2) + [.4;.5;.5;.4];
    bracketspos(1:2,2) = max(mspos(:,2)) + addonfix(1);
    bracketspos(5:6,2) = max(mspos(:,2)) + addonfix(1);
    bracketspos(3:4,2) = min(mspos(:,2)) + addonfix(3);
    bracketspos(7:8,2) = min(mspos(:,2)) + addonfix(3);
    adducttextpos = [max(bracketspos(:,1))+0.1,mean(bracketspos(:,2))];
    plotinfoout.bracketspos = bracketspos;
    plotinfoout.adducttextpos = adducttextpos;
end

end