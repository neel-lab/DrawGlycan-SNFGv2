function [attachment,cleanseq] = getseqoptions(inseq,mode)
% GETSEQOPTIONS: extract options and option values from IUPAC sequence,
% glycan or peptide.
%
% Syntax:
% attachment = getseqoptions(inseq,mode)
%
% Input:
% GLYseq: string, IUPAC condensed glycan/peptide sequence.
% mode: 'g' for glycan or 'p' for peptide.
%
% Output:
% attachment: structure, field names are option names, values are 1 x 2
% cell arrays, 1st element is option value, 2nd is which monosac. it is
% associated with
%
% Note:
% See document of GETGLYDRESSING and DRAWGLYCAN for available options.
% If there are any dashes inside quotation mark pairs, they will be
% recognized as options, therefore ruin the recognition result.
%
% Example:
% attachment = getseqoptions('Aa( -NR "x")Bb[Ff(-adduct "x34" -repeat -R)G]C[Hh(-anomer "xswe" -repeat -adduct "x25")Ii(-repeat)Jj]Dd()Ee')
%
% attachment =
%
%         NR: {'"x"'  [1]}
%          R: {[]  [3]}
%     adduct: {2x2 cell}
%     anomer: {'"xswe"'  [6]}
%     repeat: {3x2 cell}
%
% Children function:
% N/A
%

% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved

cleanseq = inseq;
switch lower(mode)
    case 'g'
        [opttype,opttypestart,opttypeend] = regexp(inseq,DrawGlycanPara.regexp_option,'match','start','end');
        [optcontent,optcontentstart,optcontentend] = regexp(inseq,DrawGlycanPara.regexp_optionvalue,'match','start','end');
        [~,MonoIndex,~,~] = Split(inseq);
        opttypemspos = zeros(size(opttypestart));
        for i = 1:length(opttypestart)
            opttypemspos(i) = find(opttypestart(i) > MonoIndex,1,'last');
        end
        allopt = [opttypestart,optcontentstart;opttypeend,optcontentend];  % matching options with contents
        if ~isempty(allopt)
            opttype = strrep(opttype,'-','');  % recognize option names by "-" prefix
            alloptstr = [opttype,optcontent];
            alloptind = [ones(size(opttypestart)),zeros(1,length(optcontent))];
            [~,sortind] = sortrows(allopt',1);
            alloptstr = alloptstr(sortind);
            alloptind = alloptind(sortind);
            whosoptind = find(alloptind);
            ind = 1;
            alloptgroup = zeros(length(opttype),2);
            for i = 1:size(alloptgroup,1)
                if whosoptind(i) < numel(alloptind)
                    if alloptind(whosoptind(i) + 1) == 0
                        alloptgroup(ind,:) = [whosoptind(i),whosoptind(i) + 1];
                    else
                        alloptgroup(ind,:) = [whosoptind(i),whosoptind(i)];
                    end
                else
                    alloptgroup(ind,:) = [whosoptind(i),whosoptind(i)];
                end
                ind = ind + 1;
            end
            if ~isempty(opttype)
                glymodinfo = [DrawGlycanPara.intglymodinfo,DrawGlycanPara.extglymodinfo,DrawGlycanPara.intglybondinfo,DrawGlycanPara.glyidentityinfo];
                attachitemind = ismember(upper(opttype),glymodinfo);
                if any(attachitemind)
                    extopttype = opttype(attachitemind);
                    extoptgroup = alloptgroup(attachitemind,:);
                    extopttypemspos = opttypemspos(attachitemind);
                    [uniattachtyp,~,attachtypeind] = unique(upper(extopttype));
                    for i = 1:max(attachtypeind)
                        optstrind = extoptgroup(attachtypeind == i,:);
                        attachcont = cell(size(optstrind,1),2);
                        tempopttypemspos = extopttypemspos(attachtypeind == i);
                        for j = 1:size(optstrind,1)
                            if diff(optstrind(j,:)) ~= 0
                                attachcont{j,1} = alloptstr{optstrind(j,2)};
                                attachcont{j,1} = attachcont{j,1}(2:end-1);
                            end
                            attachcont{j,2} = tempopttypemspos(j);
                        end
                        attachment.(upper((uniattachtyp{i}))) = attachcont;
                    end
                else
                    attachment = [];
                end
            else
                attachment = [];
            end
            for i = 1:length(optcontentstart)
                cleanseq(optcontentstart(i):optcontentend(i)) = blanks(optcontentend(i)-optcontentstart(i)+1);
            end
            [mslinkagestart,mslinkageend] = regexp(cleanseq,DrawGlycanPara.regexp_monosaclinkage,'start','end');  % based on ()
            for i = length(mslinkagestart):-1:1
                if mslinkageend(i) - mslinkagestart(i) > 1
                    templinkageseq = cleanseq(mslinkagestart(i) + 1:mslinkageend(i) - 1);
                    localopt = allopt(:,allopt(1,:) > mslinkagestart(i) & allopt(2,:) < mslinkageend(i));
                    localopt = [min(localopt(1,:)),max(localopt(2,:))] - mslinkagestart(i);
                    if ~isempty(localopt)
                        templinkageseq(localopt(1):localopt(2)) = '';
                    end
                    templinkageseq = strtrim(templinkageseq);
                    cleanseq = [cleanseq(1:mslinkagestart(i)),templinkageseq,cleanseq(mslinkageend(i):end)];
                end
            end
        else
            attachment = [];
        end
    case 'p'
        [optcontentstart,optcontentend,optcontent] = regexp(inseq,DrawGlycanPara.regexp_optionvalue,'start','end','match');  % clear " "
        for i = 1:length(optcontentstart)
            cleanseq(optcontentstart(i):optcontentend(i)) = blanks(optcontentend(i)-optcontentstart(i)+1);
        end
        [opttypestart,opttypeend,opttype] = regexp(cleanseq,DrawGlycanPara.regexp_option,'start','end','match');
        [alltypestart,alltypeend,alltype] = regexp(cleanseq,DrawGlycanPara.regexp_monosaclinkage,'start','end','match');
        for i = 1:length(alltypestart)
            cleanseq(alltypestart(i):alltypeend(i)) = blanks(alltypeend(i)-alltypestart(i)+1);
        end
        [ngmodstart,ngmodend,ngmodtype] = regexp(cleanseq,'<(.*?)>','start','end','match');
        for i = 1:length(ngmodstart)
            cleanseq(ngmodstart(i):ngmodend(i)) = blanks(ngmodend(i)-ngmodstart(i)+1);
        end
        [aaIndex,aaend,aa] = regexp(cleanseq,DrawGlycanPara.regexp_AA,'start','end','match');
        cleanseq(strfind(cleanseq,' ')) = '';  % here's the clean sequence
        opttypemspos = zeros(size(opttypestart));
        for i = 1:length(opttypestart)
            opttypemspos(i) = find(opttypestart(i) > aaIndex,1,'last');
        end
        allopt = [opttypestart,optcontentstart;opttypeend,optcontentend];  % matching options with contents
        if ~isempty(allopt)
            opttype = strrep(opttype,'-','');
            alloptstr = [opttype,optcontent];
            alloptind = [ones(size(opttypestart)),zeros(1,length(optcontent))];
            [~,sortind] = sortrows(allopt',1);
            alloptstr = alloptstr(sortind);
            alloptind = alloptind(sortind);
            whosoptind = find(alloptind);
            ind = 1;
            alloptgroup = zeros(length(opttype),2);
            for i = 1:size(alloptgroup,1)
                if whosoptind(i) < numel(alloptind)
                    if alloptind(whosoptind(i) + 1) == 0
                        alloptgroup(ind,:) = [whosoptind(i),whosoptind(i) + 1];
                    else
                        alloptgroup(ind,:) = [whosoptind(i),whosoptind(i)];
                    end
                else
                    alloptgroup(ind,:) = [whosoptind(i),whosoptind(i)];
                end
                ind = ind + 1;
            end
            if ~isempty(opttype)
                attachitemind = ismember(upper(opttype),DrawGlycanPara.aamodinfo);
                if any(attachitemind)
                    extopttype = opttype(attachitemind);
                    extoptgroup = alloptgroup(attachitemind,:);
                    extopttypemspos = opttypemspos(attachitemind);
                    [uniattachtyp,~,attachtypeind] = unique(extopttype);
                    for i = 1:max(attachtypeind)
                        optstrind = extoptgroup(attachtypeind == i,:);
                        attachcont = cell(size(optstrind,1),2);
                        tempopttypemspos = extopttypemspos(attachtypeind == i);
                        for j = 1:size(optstrind,1)
                            if diff(optstrind(j,:)) ~= 0
                                attachcont{j,1} = alloptstr{optstrind(j,2)};
                                attachcont{j,1} = attachcont{j,1}(2:end-1);
                            end
                            attachcont{j,2} = tempopttypemspos(j);
                        end
                        attachment.(uniattachtyp{i}) = attachcont;
                    end
                else
                    attachment = [];
                end
            else
                attachment = [];
            end
        else
            attachment = [];
        end
        if ~isempty(ngmodstart)
            attachcont = cell(length(ngmodstart),2);
            for i = 1:length(ngmodstart)
                ngmodmspos = find(ngmodstart(i) > aaIndex,1,'last');
                attachcont(i,:) = {ngmodtype{i}(2:end-1),ngmodmspos};
            end
            attachment.NGMOD = attachcont;
        end
end
end