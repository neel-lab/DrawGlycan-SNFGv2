# DrawGlycan-SNFGv2
Open source code for DrawGlycan-SNFG (version 2)

This is a fully functional version of DrawGlycan-SNFG (version 2). To use it:
1. Copy these files to a folder
2. Open MATLAB and navigate to the folder
3. Then draw the glycan using a command like:<br/>
to draw Glucose<br/>
     \>> drawglycan('Glc')            <br/>
to draw a LacNAc chain *etc*.<br/>
    \>> drawglycan('Gal(b4)GlcNAc(b?)') 
to draw multiple glycans in a single window<br/>
    \>> drawglycanTile({'Gal(b3)GalNAc','Gal(b3)GlcNAc','Gal','Gal(b4)GlcNAc'})<br/>
  There are many more examples of IUPAC strings you can use at http://virtualglycome.org/DrawGlycan/DrawExamples.html
 4. To deploy DrawGlycan-SNFG in PYTHON environment, use the matlab.engine API: https://www.mathworks.com/help/matlab/matlab-engine-for-python.html
 5. To get help, please see the documentation at http://virtualglycome.org/DrawGlycan/DrawGlycanSNFG2_Manual.pdf
   or simply write <br/>
    \>> help drawglycan
