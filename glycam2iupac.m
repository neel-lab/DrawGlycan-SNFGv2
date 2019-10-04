function IUPACout = glycam2iupac(glycam)
% GLYCAM2IUPAC: translate glycan sequence from GLYCAM to IUPAC condensed
%
% Syntax:
% IUPACout = glycam2iupac(glycam)
%
% Input:
% glycan: glycan sequence in GLYCAM format.
%
% Output:
% IUPACout: glycan sequence in IUPAC format.
%
% Note:
% UNAVAILABLE - STILL IN DEVELOPMENT
%
% Example:
% N/A
%
% Children function:
% N/A

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


IUPACout = glycam;
bondlinks = unique(regexp(glycam,'(a|b)[12]-[0-9|OH|OME|OtBu]+','match'));
for i = 1:length(bondlinks)
    bondlinksnew = ['(',bondlinks{i},')'];
    IUPACout = strrep(IUPACout,bondlinks{i},bondlinksnew);
end
specmsexp = 'DNeup5Ac|DNeup5Gc|KDN|KDO|LGalpNAc|DGalpNAc|LGlcpNAc|DGlcpNAc|LManpNAc|DManpNAc|LGalpA|DGalpA|LGlcpA|DGlcpA|LIdopA|DIdopA';
specms = unique(regexp(glycam,specmsexp,'match'));
specmsori = {'DNeup5Ac','DNeup5Gc','KDN','KDO','LGalpNAc','DGalpNAc','LGlcpNAc','DGlcpNAc','LManpNAc','DManpNAc','LGalpA','DGalpA','LGlcpA','DGlcpA','LIdopA','DIdopA'};
specmsnew = {'Neu5Ac','Neu5Gc','KDN','KDO','GalNAc','GalNAc','GlcNAc','GlcNAc','ManNAc','ManNAc','GalA','GalA','GlcA','GlcA','IdoA','IdoA'};
specms = specmsori(ismember(specmsori,specms));
newspecms = specmsnew(ismember(specmsori,specms));
for i = 1:length(specms)
    IUPACout = strrep(IUPACout,specms{i},newspecms{i});
end
regmsexp = 'Man|Gal|Glc|Ido|All|Alt|Gul|Tal|Xyl|Lyx|Rib|Ara|Fru|Psi|Sor|Tag|Fuc|Rha|Qui';
regms = unique(regexp(glycam,['(L|D)(',regmsexp,')(f|p)'],'match'));
for i = 1:length(regms)
    IUPACout = strrep(IUPACout,regms{i},regms{i}(2:end-1));
end
end