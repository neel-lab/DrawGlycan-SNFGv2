function drawglycanTile(glycanlist,varargin)
% INPUT: 
% glycanlist: cell array of strings, each containing glycan sequence
% optional inputs:
% 1. rownum: number of rows for each subplot
% 2. colnum: number of columns for each subplot
% example: drawglycanTile({'Gal(b3)GalNAc','Gal(b3)GlcNAc','Gal','Gal(b4)GlcNAc'}) 

if nargin>1
    rownum=varargin{1};
    colnum=varargin{2};
else
    rownum=ceil(length(glycanlist)^0.5);
    colnum=rownum;
end
glycanlist = glycanlist(:);
glycanlist = glycanlist(~cellfun(@isempty,glycanlist));
for ii = 1:rownum*colnum:length(glycanlist)
    tempglycanlist = glycanlist(ii:min(ii + rownum * colnum - 1,length(glycanlist)));
    figure('Name',['Glycan No. ',num2str(ii),' - ',num2str(min(ii + rownum * colnum - 1,length(glycanlist)))],'NumberTitle','off');
    for jj = 1:length(tempglycanlist)
        ax = subplot(rownum,colnum,jj);
        hold(ax,'on');
        drawglycan(tempglycanlist{jj},'inputformat','IUPAC','figurehandle',ax);
        set(ax,'visible',false);
    end
end
end