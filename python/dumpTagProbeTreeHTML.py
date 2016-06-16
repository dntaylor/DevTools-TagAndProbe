#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse, json, pickle, os, re, subprocess

from DevTools.Utilities.utilities import python_mkdir

__gitversion__ = subprocess.check_output(["git", "describe", "--always"]).strip()

def rootFileType(string) :
    file = ROOT.TFile.Open(string)
    if not file :
        raise argparse.ArgumentTypeError(string+' could not be opened!')
    return file

def get2DPlot(effDir) :
    canvas = None
    plotsDir = effDir.Get('fit_eff_plots')
    if not plotsDir :
        return
    for key in plotsDir.GetListOfKeys() :
        if re.match('probe_.*eta.*_probe_.*_PLOT.*', key.GetName()) :
            canvas = key.ReadObj()
            break

    if not canvas :
        raise Exception('No canvas found in %s' % effDir.GetName())
    
    # Stylize?
    canvas.SetLogy(True)
    plot = canvas.GetPrimitive(canvas.GetName())

    ROOT.SetOwnership(canvas, False)
    return canvas

def parseFitTree(baseDirectory, outputDirectory) :
    '''
        Generates HTML files in outputDirectory
    '''

    python_mkdir(outputDirectory)
    with open(os.path.join(outputDirectory, 'index.html'), 'w') as index :
        index.write(HTML.header.format(
            rootFile=baseDirectory.GetFile().GetName(),
            baseDirectory=baseDirectory.GetName()
        ))

        for effDir in subDirs(baseDirectory) :
            effName = effDir.GetName()

            effoutDir = os.path.join(outputDirectory, effName)
            python_mkdir(effoutDir)

            parseEfficiencyDir(effDir, effoutDir, index)
        
        index.write(HTML.footer.format(
                version=__gitversion__
            ))

def parseEfficiencyDir(effDir, outputDirectory, index) :
    effName = effDir.GetName()
    rows = ''
    codeRows = []

    for effBinDir in fitSubDirs(effDir) :
        binName = effBinDir.GetName()

        if not effBinDir.Get('fitresults') :
            continue

        (row, code) = parseEfficiencyBin(effBinDir, outputDirectory)
        codeRows.append(code)
        rows += row

    plot = get2DPlot(effDir)
    plotName = 'efficiency2DPtEta.png'
    if plot :
        plot.Print(os.path.join(outputDirectory, plotName))

    macroVariables = effDir.Get('variables').GetTitle()
    with open(os.path.join(outputDirectory, effName+'.C'), 'w') as fout :
        fout.write("// Made with DevTools/TagAndProbe version: %s\n" % __gitversion__)
        fout.write("""
enum Variation {
    CENTRAL,
    STAT_UP,
    STAT_DOWN,
    SYST_ALT_TEMPL,
    SYST_TAG_PT30,
    SYST_CMSSHAPE,
    EFF_DATA,
    EFF_DATA_ERRSYM,
    EFF_MC,
    EFF_MC_ERRSYM,
};
""")
        fout.write("float %s(%s, Variation variation) {\n" % (effName, macroVariables))
        for row in codeRows :
            fout.write(row+"\n")
        fout.write("    return 1.; // Default\n")
        fout.write("}\n")

    index.write(HTML.effDirItem.format(
        effName=effName,
        eff2DPtEtaImage=os.path.join(effName, plotName),
        tableRows=rows
    ))

def makeChi2(rooPad, nFitParam) :
    # This is so fucking stupid
    data = rooPad.FindObject('h_data_binned')
    model = rooPad.FindObject('pdfPass_Norm[mass]')
    chi2 = 0.
    for i in range(data.GetN()) :
        x = ROOT.Double()
        y = ROOT.Double()
        data.GetPoint(i, x, y)
        yerr = data.GetErrorY(i)
        ypred = model.interpolate(x)
        if y > 0 :
            chi2 += (y-ypred)**2 / yerr**2
    return (chi2, data.GetN()-nFitParam)

def makeChi2_roofit(rooPad, nFitParam) :
    data = rooPad.FindObject('h_data_binned')
    model = rooPad.FindObject('pdfPass_Norm[mass]')
    Redchi2 = model.chiSquare(data, nFitParam)
    dof = data.GetN() - nFitParam
    return (Redchi2*dof, dof)

def parseEfficiencyBin(effBinDir, outputDirectory) :
    fitResults = effBinDir.Get('fitresults')
    # https://root.cern.ch/root/html/tutorials/roofit/rf607_fitresult.C.html
    effValue = fitResults.floatParsFinal().find('efficiency')
    dataEff = effValue.getVal()
    dataEffErrHi = effValue.getErrorHi()
    dataEffErrLo = effValue.getErrorLo()

    if effValue.getVal()+effValue.getErrorHi() > 1. :
        print "Found one! :"
        effValue.Print()

    canvas = effBinDir.Get('fit_canvas')
    passing = canvas.FindObject('fit_canvas_1')
    firstPlot = passing.GetListOfPrimitives()[0]
    firstPlot.SetTitle('Passing '+effBinDir.GetName())
    fullPlot = os.path.join(outputDirectory, effBinDir.GetName()+'.png')
    canvas.Print(fullPlot)
    smallPlot = fullPlot.replace('.png','_small.png')
    subprocess.call(['convert', '-gravity', 'north', '-crop', '100%x50%', fullPlot, smallPlot])

    cHist = fitResults.correlationHist()
    cHist.SetTitle('')
    cHist.SetContour(256)
    canvasCorr = ROOT.TCanvas('correlation', '', 800,600)
    canvasCorr.SetMargin(.2,.12,.1,.01)
    cHist.Draw('colz')
    canvasCorr.Print(fullPlot.replace('.png','_correlation.png'))
    cHist.SetDirectory(0)
    del cHist

    nll = fitResults.minNll()

    mcHist = effBinDir.Get('mc_cutCount')
    mcPass = mcHist.GetBinContent(1)
    mcTotal= mcHist.GetBinContent(2)
    mcEff = mcPass/mcTotal if mcTotal else 0.
    mcEffLo = ROOT.TEfficiency.ClopperPearson(int(mcTotal), int(mcPass), 0.68, False)
    mcEffHi = ROOT.TEfficiency.ClopperPearson(int(mcTotal), int(mcPass), 0.68, True)

    scaleFactor = dataEff / mcEff if mcEff else 0.
    maxSf = (dataEff+dataEffErrHi)/mcEffLo if mcEffLo else 0.
    scaleFactorErrHi = maxSf-scaleFactor
    minSf = (dataEff+dataEffErrLo)/mcEffHi if mcEffHi else 0.
    scaleFactorErrLo = minSf-scaleFactor

    with open(os.path.join(outputDirectory, effBinDir.GetName()+'.parameters.txt'), 'w') as paramOut :
        #paramOut.write("# Fit chi2/dof: %f / %d\n" % (chi2, dof))
        params = fitResults.floatParsFinal()
        for p in xrange(params.getSize()):
            myPar = params.at(p)
            paramOut.write("%s[%.3f,%.3f,%.3f]\n"%(myPar.GetName(), myPar.getVal(), myPar.getMin(), myPar.getMax()))

    colors = ['#79ff00', '#ffff00', '#ff7700', '#ff0000']
    thresholds = [1., 5., 10., 100.]
    fitStatusColor = '#00ff00' 
    for i in range(4) :
        if nll > thresholds[i] :
            fitStatusColor = colors[i]

    binName = effBinDir.GetName()

    row = HTML.effTableRow.format(
            fitStatusColor=fitStatusColor,
            efficiencyNice='%.4f +%.4f %.4f' % (dataEff, dataEffErrHi, dataEffErrLo),
            mcEfficiency='%.4f (%d/%d)' % (mcEff, mcPass, mcTotal),
            scaleFactor='%.4f +%.4f %.4f' % (scaleFactor, scaleFactorErrHi, scaleFactorErrLo),
            effName=effBinDir.GetDirectory('..').GetName(),
            binName=binName,
            binNameNice=re.sub('__pdfSignal.*', '', binName),
            numSignalAll=fitResults.floatParsFinal().find('numSignalAll').getVal(),
            nll=nll,
            fitStatus=':'.join(['%d' % fitResults.statusCodeHistory(i) for i in range(fitResults.numStatusHistory())]),
        )
    row += '''<tr><td colspan="7" style="text-align: center;"><img src="{effName}/{binName}_small.png"></td></tr>'''.format(
            effName=effBinDir.GetDirectory('..').GetName(),
            binName=binName,
        )

    codeRow = effBinDir.Get('cutString').GetTitle()

    return (row, codeRow)

def subDirs(dir) :
    for key in dir.GetListOfKeys() :
        if key.IsFolder() :
            yield key.ReadObj()

def fitSubDirs(dir) :
    for key in dir.GetListOfKeys() :
        if key.IsFolder() :
            yield key.ReadObj()

# HTML Templates
class HTML :
    header = '''<html>
<head>
    <title>Efficiency summary for {rootFile}/{baseDirectory}</title>
    <style type="text/css">
        div.effResult {{
            margin: auto;
            width: 80%;
            border: 2px solid #888888;
        }}

        table {{
            width: 100%;
        }}
    </style>
</head>
<body>
<h2>Efficiency summary for {rootFile}/{baseDirectory}</h2>
'''

    effDirItem = '''<div class="effResult">
    <h3>Efficiency results for {effName}</h3><br />
    Function: <a href="{effName}/{effName}.C">{effName}.C</a><br />
    Latex: <a href="{effName}/table.tex">table.tex</a><br />
    <div style="text-align: center;"><b>Summary Plot</b><br /><img src="{effName}/scaleFactor_vs_pt.png" /><br /><img src="{effName}/scaleFactor_vs_eta.png" /></div>
    <table>
        <tr style="background: #cccccc">
            <td>Bin name</td>
            <td>Data Efficiency</td>
            <td>MC Cut&Count</td>
            <td>Scale Factor</td>
            <td>Fit status <a href="https://root.cern.ch/root/htmldoc/ROOT__Minuit2__Minuit2Minimizer.html#ROOT__Minuit2__Minuit2Minimizer:Minimize">?</a></td>
            <td>Fit signal yield</td>
            <td>Extras</td>
        </tr>
        {tableRows}
    </table>
</div>'''

    effTableRow = '''
        <tr style="background: {fitStatusColor}">
            <td><a href="{effName}/{binName}.png">{binNameNice}</a></td>
            <td>{efficiencyNice}</td>
            <td>{mcEfficiency}</td>
            <td>{scaleFactor}</td>
            <td>{fitStatus}</td>
            <td>{numSignalAll:.0f}</td>
            <td><a href="{effName}/{binName}.parameters.txt">Params</a> <a href="{effName}/{binName}_correlation.png">Corr.</a> NLL: {nll:.0f}</td>
        </tr>'''


    footer = '''</table>
Generated using DevTools/TagAndProbe version {version}<br />
</body>
</html>
'''

def setPalette() :
    import array
    def readColors(fileName) :
        with open(fileName) as palette :
            for i, line in enumerate(palette) :
                (r,g,b) = map(float, line.strip().split(','))
                c = ROOT.TColor(1000+i,r,g,b)
                ROOT.SetOwnership(c, False)
                yield c.GetNumber()

    # http://peterkovesi.com/projects/colourmaps/ColourmapTheory/index.html
    colors = array.array('i', readColors('{0}/src/DevTools/TagAndProbe/data/diverging_bwr_55-98_c37_n256.csv'.format(os.environ['CMSSW_BASE'])))
    ROOT.TColor.SetPalette(len(colors), colors)
    ROOT.gStyle.SetNumberContours(len(colors))

def main() :
    parser = argparse.ArgumentParser(description='Dumps fit info generated by TagProbeFitTreeAnalyzer into HTML summary')
    parser.add_argument('--data', help='Data fit tree name', type=rootFileType)
    parser.add_argument('--output', '-o', help='Directory name for output', required=True)
    parser.add_argument('--input', '-i', help='Directory name in root files to load', required=True)
    args = parser.parse_args()
    
    ROOT.gStyle.SetPaintTextFormat('0.4f')
    ROOT.gStyle.SetOptTitle(True)
    setPalette()

    python_mkdir(args.output)

    parseFitTree(args.data.Get(args.input), args.output)

if __name__ == '__main__' :
    main()

