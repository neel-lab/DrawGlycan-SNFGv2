classdef MonoColor
    %% MONOCOLOR: color definition for monosaccharides
    % Note: 6dAlt and 6dTal are presented as dAlt6 and dTal6 to avoid variable
    % naming issues.
    
    %
    % DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
    %(c) 2019, Research Foundation for State University of New York. All rights reserved
    %
    
    
    properties (Constant)
        % Color RGB values
        White = [1 1 1];
        Blue = [0 .565 .737];
        Green = [0 .651 .318];
        Yellow = [1 .831 0];
        LightBlue = [.561 .8 .914];
        Pink = [0.965 0.620 0.631];
        Purple = [.647 .263 0.6];
        Brown = [.631 .478 .302];
        Orange = [.957 .475 .125];
        Red = [.929 .11 .141];
        Black = [0 0 0];
        None = [0 0 0];
        
        % Hexose
        Hex = 'White';
        Hexose = 'White';
        Glc = 'Blue';
        Man = 'Green';
        Gal = 'Yellow';
        Gul = 'Orange';
        Alt = 'Pink';
        All = 'Purple';
        Tal = 'LightBlue';
        Ido = 'Brown';
        
        % HexNAc
        HexNAc = 'White';
        GlcNAc = 'Blue';
        ManNAc = 'Green';
        GalNAc = 'Yellow';
        GulNAc = 'Orange';
        AltNAc = 'Pink';
        AllNAc = 'Purple';
        TalNAc = 'LightBlue';
        IdoNAc = 'Brown';
        
        % HexN
        Hexosamine = 'White';
        GlcN = 'Blue';
        ManN = 'Green';
        GalN = 'Yellow';
        GulN = 'Orange';
        AltN = 'Pink';
        AllN = 'Purple';
        TalN = 'LightBlue';
        IdoN = 'Brown';
        
        % HexA
        Hexuronate = 'White';
        GlcA = 'Blue';
        ManA = 'Green';
        GalA = 'Yellow';
        GulA = 'Orange';
        AltA = 'Pink';
        AllA = 'Purple';
        TalA = 'LightBlue';
        IdoA = 'Brown';
        
        % Deoxyhexose
        Deoxyhexose = 'White';
        Qui   = 'Blue';
        Rha   = 'Green';
        Fuc   = 'Red';
        
        % DeoxyhexNAc
        DeoxyhexNAc = 'White';
        QuiNAc   = 'Blue';
        RhaNAc   = 'Green';
        FucNAc   = 'Red';
        
        % Di-Deoxyhexose
        Dideoxyhexose = 'White';
        Oli   = 'Blue';
        Tyv   = 'Green';
        Abe   = 'Orange';
        Par   = 'Pink';
        Dig   = 'Purple';
        Col   = 'LightBlue';
        
        % Pentose
        Pentose = 'White';
        Ara   = 'Green';
        Lyx   = 'Yellow';
        Xyl   = 'Orange';
        Rib   = 'Pink';
        
        % Nonulosonate
        Deoxynonulosonate = 'White';
        Kdn    = 'Green';
        Neu5Ac = 'Purple';
        Neu5Gc = 'LightBlue';
        Neu    = 'Brown';
        Sia = 'Red';
        
        % Dideoxynonulosonate
        Dideoxynonulosonate = 'White'
        Pse = 'Green';
        Leg = 'Yellow';
        Aci = 'Pink';
        
        % Unknown
        Unknown = 'White';
        Bac      = 'Blue';
        LDManHep = 'Green';
        Kdo      = 'Yellow';
        Dha      = 'Orange';
        DDManHep = 'Pink';
        MurNAc   = 'Purple';
        MurNGc   = 'LightBlue';
        Mur      = 'Brown';
        
        % Assigned
        Assigned = 'White';
        Api = 'Blue';
        Fru = 'Green';
        Tag = 'Yellow';
        Sor = 'Orange';
        Psi = 'Pink';
        NonGly = 'Black';
        Blank = 'None';
    end
end