import ROOT
import json
from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from variables.drs import calibrateDRSPeakTS
from utils.html_generator import generate_html
from utils.parser import get_args
from plotting.my_function import DrawHistos, LHistos2Hist
from utils.timing import auto_timer
from utils.root_setup import setup_root
from utils.utils import number_to_string
import numpy as np
import matplotlib.pyplot as plt
import re
import os

setup_root(n_threads=10, batch_mode=True, load_functions=True)
#input_file = "/lustre/work/jweijie/photonstudy/run1635_260424121214_converted.root"
#input_file = "/lustre/work/jweijie/photonstudy/run1636_260424134933_converted.root"
#input_file = "/lustre/work/jweijie/photonstudy/run1667_260501172036_converted.root"
input_file = "/lustre/work/jweijie/photonstudy/run1672_260501191124_converted.root"
#input_file = "/lustre/work/jweijie/photonstudy/run1666_260501155847_converted.root"
run_number = 1672
path = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results"
rootdir = f"{path}/root/Run{run_number}/"
plotdir = f"{path}/plots/Run{run_number}/"
os.makedirs(rootdir, exist_ok=True)
os.makedirs(plotdir, exist_ok=True)

#DRS_ana = [
#    "DRS_Board0_Group2_Channel8",
#    "DRS_Board0_Group2_Channel2"
#
#]
#DRS_channels = [
#    "DRS_Board0_Group2_Channel2"
#]
DRS_ana = [
    "DRS_Board2_Group0_Channel8",
    "DRS_Board2_Group0_Channel1"

]
DRS_channels = [
    "DRS_Board2_Group0_Channel1"
]
#DRS_ana = []
#DRS_channels = []
#for b in range(4):
#    for g in range(4):
#        for c in range(9):
#            DRS_ana.append(f"DRS_Board{b}_Group{g}_Channel{c}")
#            if c != 8 :
#                DRS_channels.append(f"DRS_Board{b}_Group{g}_Channel{c}")



rdf = ROOT.RDataFrame("EventTree", input_file)

def processChannel(rdf, drs_ana, drs_chan):
    hists = []
    rdf = rdf.Define("TS", "FillIndices(1024)")
    for varname in drs_ana:
        rdf = rdf.Define(
            f"{varname}_bl",
            f"compute_median({varname})"
        )
        rdf = rdf.Define(
                f"{varname}_blsub",
                f"{varname} - {varname}_bl"
            )
        rdf = rdf.Define(f"{varname}_PeakTS",
                             f"ArgMinRange({varname}_blsub, 100, 500, -600.0)")


    for chan in drs_chan:
        channel_TS = re.sub(r"_Channel[0-7]", "_Channel8", chan)
        hist2d_DRS_VS_TS_trigger = rdf.Histo2D((
                    f"hist2d_DRS_vs_TS_{channel_TS}_raw", "", 1024, 0, 1024, 100, -1000, 1000), "TS", f"{channel_TS}_blsub")
        
        rdf = rdf.Define(
                f"{chan}_AlignedTS", f"TS - (int){channel_TS}_PeakTS"
            )

        #rdf = rdf.Define(
        #        f"{chan}_Sum", f"SumRange({chan}_blsub, 660, 780)")
        #rdf = rdf.Define(
        #        f"{chan}_Sum", f"SumRange({chan}_blsub, 0, 1024)")
        rdf = rdf.Define(
                f"{chan}_Sum", f"SumRange({chan}_blsub, 0, 1024)")
        h_DRS_sum= rdf.Histo1D((
                    f"hist_DRS_sum_{chan}",
                    f"DRS sum;#ADC;Counts",
                    500, -1e3, 7e3),
                    f"{chan}_Sum"
                )
        rdf = rdf.Define(
                f"{chan}_blsub_square",
                f"{chan}_blsub*{chan}_blsub"
            )
        rdf = rdf.Define(
                f"{chan}_RMS",
                f"sqrt(SumRange({chan}_blsub_square, 0, 1024)/1024)"
            )
        h_DRS_RMS= rdf.Histo1D((
                    f"hist_DRS_RMS_{chan}",
                    f"DRS RMS;#ADC;Counts",
                    100, 0, 100),
                    f"{chan}_RMS"
                )
        hist2d_DRS_VS_TS_raw = rdf.Histo2D((
                    f"hist2d_DRS_vs_TS_{chan}_raw", "", 1024, 0, 1024, 100, -1000, 1200), "TS", f"{chan}_blsub")
        hist2d_DRS_VS_TS = rdf.Histo2D((
                    f"hist2d_DRS_vs_TS_{chan}", "", 1000, -100, 900, 140, -200, 1200), f"{chan}_AlignedTS", f"{chan}_blsub")
        hists.append(hist2d_DRS_VS_TS)
        hists.append(hist2d_DRS_VS_TS_raw)
        hists.append(hist2d_DRS_VS_TS_trigger)
        hists.append(h_DRS_sum)
        hists.append(h_DRS_RMS)
        #rdf_filtered = rdf.Filter(
        #    f"{chan}_Sum>2650"
        #)
        #for ievent in range(10):
        #    hist2d_DRS_VS_TS_entry = rdf_filtered.Range(ievent,ievent+1).Histo2D((
        #            f"hist2d_DRS_vs_TS_{chan}_entry{ievent}", "", 1024, 0, 1024, 100, -20, 600), "TS", f"{chan}_blsub")
        #    hists.append(hist2d_DRS_VS_TS_entry)

    return hists
def makeDRSPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRS"
    infile_name = f"{rootdir}/drs.root"
    infile = ROOT.TFile(infile_name, "READ")
    for hname in [f"hist2d_DRS_vs_TS_{varname}" for varname in DRS_channels] + [f"hist2d_DRS_vs_TS_{varname}_raw" for varname in DRS_ana]:
        h = infile.Get(hname)
        out_name = hname.replace("hist2d_","")
        xtitle = h.GetXaxis().GetTitle()
        ytitle = h.GetYaxis().GetTitle()
        xtitle = "TS"
        ytitle = "#ADC"
        if "raw" in hname:
            xmin = 0
            xmax = 1024
            ymin = -1000
            ymax = 800
        else:
            xmin = -100
            xmax = 900
            ymin = -200
            ymax = 800
        DrawHistos([h], "", xmin, xmax, xtitle, ymin, ymax, ytitle,
                    out_name,
                    dology=False, addOverflow=False, addUnderflow=False, doth2=True, drawoptions="COLZ", zmin=1, zmax=1e4, dologz=True,
                    outdir=outdir_plots, run_number=run_number)
    '''
    for hname in [f"hist_DRS_sum_{varname}" for varname in DRS_channels]:
        h = infile.Get(hname)
        out_name = hname.replace("hist_","")
        xtitle = h.GetXaxis().GetTitle()
        ytitle = h.GetYaxis().GetTitle()
        xmin = -1e3
        xmax = 7e3
        DrawHistos([h], ["ADC sum"], xmin, xmax, xtitle, 1, None, ytitle,
            out_name,
            dology=True, drawoptions="HIST", mycolors=[2], addOverflow=False, addUnderflow=False,
            outdir=outdir_plots, run_number=run_number)
    
    for hname in [f"hist_DRS_RMS_{varname}" for varname in DRS_channels]:
        h = infile.Get(hname)
        out_name = hname.replace("hist_","")
        xtitle = h.GetXaxis().GetTitle()
        ytitle = h.GetYaxis().GetTitle()
        xmin = 0
        xmax = 100
        DrawHistos([h], ["RMS"], xmin, xmax, xtitle, 1, None, ytitle,
            out_name,
            dology=True, drawoptions="HIST", mycolors=[2], addOverflow=False, addUnderflow=False,
            outdir=outdir_plots, run_number=run_number)
    
    '''
    for hname in [f"hist_DRS_sum_{varname}" for varname in DRS_channels]:
        h = infile.Get(hname)
        out_name = hname.replace("hist_","")
        xtitle = h.GetXaxis().GetTitle()
        ytitle = h.GetYaxis().GetTitle()
        xmin = -1e3
        xmax = 7e3
        #fit_center = [0, 500, 1000, 1450, 1900, 2400, 2900]
        #fit_width = [100, 100, 100, 100, 100, 100, 100]
        #fit_amp = [250, 400, 400, 300, 200, 150, 60]
        #fit_center = [0, 1000, 1800, 2700, 3400, 4000, 4800]
        #fit_width = [300, 300, 300, 300, 300, 300, 300]
        #fit_amp = [70, 150, 150, 150, 100, 80, 60]
        #fit_center = [100, 850, 1800, 2700, 3500, 4400, 5100]
        #fit_width = [300, 300, 300, 300, 300, 300, 300]
        #fit_amp = [600, 1000, 800, 600, 400, 200, 100]
        #fit_color = [3,4,6,7,8,9,46]
        fit_center = [0, 900, 1800, 2800, 3600, 4400]
        fit_width = [400, 400, 400, 400, 400, 400]
        fit_amp = [100, 300, 400, 350, 300, 200]
        fit_color = [3,4,6,7,8,9]
        #fitmin = -500
        #fitmax = 3150
        #fitmin = -400
        #fitmax = 5300
        #fitmin = -500
        #fitmax = 5500
        fitmin = -300
        fitmax = 5000

        #fit_func = ROOT.TF1("fit", "gaus(0) + gaus(3) + gaus(6) + gaus(9) + gaus(12) + gaus(15) + gaus(18)", fitmin, fitmax)
        fit_func = ROOT.TF1("fit", "gaus(0) + gaus(3) + gaus(6) + gaus(9) + gaus(12) + gaus(15)", fitmin, fitmax)
        fit_func.SetLineColor(1)
        
        for i,(center, width, amp) in enumerate(zip(fit_center, fit_width, fit_amp)):
            fit_func.SetParameter(i*3, amp)
            fit_func.SetParameter(i*3+1, center)
            fit_func.SetParameter(i*3+2, width)
            fit_func.SetParLimits(i*3,amp*0.1,amp*2)
            fit_func.SetParLimits(i*3+1,center-250, center+250)
            fit_func.SetParLimits(i*3+2,1,width+180)
        h.Fit(fit_func, "S R")
        a = [
            fit_func.GetParameter(0),
            fit_func.GetParameter(3),
            fit_func.GetParameter(6),
            fit_func.GetParameter(9),
            fit_func.GetParameter(12),
            fit_func.GetParameter(15),
            #fit_func.GetParameter(18),
            ]
        mu = [
            fit_func.GetParameter(1),
            fit_func.GetParameter(4),
            fit_func.GetParameter(7),
            fit_func.GetParameter(10),
            fit_func.GetParameter(13),
            fit_func.GetParameter(16),
            #fit_func.GetParameter(19)
        ]
        sigma = [
            fit_func.GetParameter(2),
            fit_func.GetParameter(5),
            fit_func.GetParameter(8),
            fit_func.GetParameter(11),
            fit_func.GetParameter(14),
            fit_func.GetParameter(17),
            #fit_func.GetParameter(20)
        ]
        fit_individuals = [
            ROOT.TF1("g1","gaus",fitmin,fitmax),
            ROOT.TF1("g2","gaus",fitmin,fitmax),
            ROOT.TF1("g3","gaus",fitmin,fitmax),
            ROOT.TF1("g4","gaus",fitmin,fitmax),
            ROOT.TF1("g5","gaus",fitmin,fitmax),
            ROOT.TF1("g6","gaus",fitmin,fitmax),
            #ROOT.TF1("g7","gaus",fitmin,fitmax)
        ]
        for g,amp,center,width,c in zip(fit_individuals,a,mu,sigma,fit_color):
            g.SetParameter(0,amp)
            g.SetParameter(1,center)
            g.SetParameter(2,width)
            g.SetLineColor(c)

        extraToDraw = ROOT.TPaveText(0.55, 0.45, 0.9, 0.75, "NDC")
        #extraToDraw = ROOT.TPaveText(0.15, 0.7, 0.7, 0.9, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.03)
        #extraToDraw.AddText(f"amp = {a[0]:.0f},{a[1]:.0f},{a[2]:.0f},")
        #extraToDraw.AddText(f"{a[3]:.0f},{a[4]:.0f},{a[5]:.0f},{a[6]:.0f}")
        #extraToDraw.AddText(f"mu = {mu[0]:.0f},{mu[1]:.0f},{mu[2]:.0f}, ")
        #extraToDraw.AddText(f"{mu[3]:.0f},{mu[4]:.0f},{mu[5]:.0f},{mu[6]:.0f}")
        #extraToDraw.AddText(f"sigma = {sigma[0]:.1f},{sigma[1]:.1f},{sigma[2]:.1f},")
        #extraToDraw.AddText(f"{sigma[3]:.1f},{sigma[4]:.1f},{sigma[5]:.1f},{sigma[6]:.1f}")
        extraToDraw.AddText(f"amp = {a[0]:.0f},{a[1]:.0f},{a[2]:.0f},")
        extraToDraw.AddText(f"{a[3]:.0f},{a[4]:.0f},{a[5]:.0f}")
        extraToDraw.AddText(f"mu = {mu[0]:.0f},{mu[1]:.0f},{mu[2]:.0f}, ")
        extraToDraw.AddText(f"{mu[3]:.0f},{mu[4]:.0f},{mu[5]:.0f}")
        extraToDraw.AddText(f"sigma = {sigma[0]:.1f},{sigma[1]:.1f},{sigma[2]:.1f},")
        extraToDraw.AddText(f"{sigma[3]:.1f},{sigma[4]:.1f},{sigma[5]:.1f}")
        DrawHistos([h], ["ADC sum"], xmin, xmax, xtitle, 1, None, ytitle,
            out_name,
            dology=True, drawoptions="HIST", mycolors=[2], addOverflow=False, addUnderflow=False, extraToDraw=fit_individuals+[fit_func,extraToDraw],
            outdir=outdir_plots, run_number=run_number)

        #nph = range(7)
        nph = range(6)
        m, b = np.polyfit(nph, mu, 1)
        plt.scatter(nph, mu, color='blue', label='mu')
        plt.plot(nph, m*np.array(nph) + b, color='red', label=f'Linear Fit: #ADC={m:.2f}nph+{b:.2f}')
        plt.xlabel("nph")
        plt.ylabel('#ADC')
        plt.title('Linearity check of #ADC versus #photons')
        plt.legend()
        plt.savefig(f"{outdir_plots}/fit_ADC_nph.png", dpi=300, bbox_inches='tight')
        

def main():
    hists = processChannel(rdf,DRS_ana, DRS_channels)
    outfile_DRS = ROOT.TFile(
            f"{rootdir}/drs.root", "RECREATE")

    for h in hists:
        h.SetDirectory(outfile_DRS)
        h.Write()

    outfile_DRS.Close()

    makeDRSPlots()
if __name__ == "__main__":
    main()


        
