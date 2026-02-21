import ROOT
import json
from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from configs.plot_config import getRangesForFERSEnergySums
from variables.drs import calibrateDRSPeakTS, getDRSPeak
from utils.html_generator import generate_html
from utils.parser import get_args
from plotting.my_function import DrawHistos, LHistos2Hist
from utils.timing import auto_timer
from utils.root_setup import setup_root
from utils.utils import number_to_string

run_number = 1513
benergy = 40
run_number_text = "1501,1507,1511,1513"
plotdir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/plots/RunElectron90"
rootdir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/root/RunElectron90"
htmldir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/html/RunElectron90"
GainCalibs = [("HG", False), ("LG", False), ("Mix", True)]
HE = (benergy >= 50)

TSmin = -80
TSmax = -20
TSCermin = -80
TSCermax = -50


def makeDRSPeakTSVSEnergyPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeakTS_VS_Energy"
    infile_name = f"{rootdir}/drspeak_vs_fers.root"
    infile = ROOT.TFile(infile_name, "READ")
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)
        for cat in ["cer", "sci"]:
            hist_names = [
                f"hist_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}",
                f"hist_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}",
                f"hist_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}",
                ]
            xtitles = [
                "Cer Quartz peak TS",
                "Cer Plastic peak TS",
                "Sci peak TS"
                ]
            ytitle = f"FERS Energy {cat} {gain}"
            for hist_name, xtitle in zip(hist_names, xtitles):
                output_name = hist_name.replace("hist_","")
                hist = infile.Get(hist_name)
                extraToDraw = ROOT.TPaveText(0.20, 0.85, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(f"correlation = {hist.GetCorrelationFactor():.3f}")
                DrawHistos([hist], "", TSmin if "Sci" in xtitle else TSCermin, TSmax if "Sci" in xtitle else TSCermax, xtitle, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"], ytitle,
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
                       outdir=outdir_plots, addOverflow=False, run_number=run_number_text,extraToDraw=[extraToDraw])
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/FERSvsDRS//EnergySum_VS_DRSTS.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html
    
def makeDRSPeakTSVSCSratioPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeakTS_VS_CSratio"
    infile_name = f"{rootdir}/drspeak_vs_fers_csratio.root"
    infile = ROOT.TFile(infile_name, "READ")
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)
        for cat in ["cer", "sci"]:
            hist_names = [
                f"hist_DRSPeakTS_Cer_Quartz_VS_FERS_CSratio_{gain}",
                f"hist_DRSPeakTS_Cer_Plastic_VS_FERS_CSratio_{gain}",
                f"hist_DRSPeakTS_Sci_VS_FERS_CSratio_{gain}",
                ]
            xtitles = [
                "Cer Quartz peak TS",
                "Cer Plastic peak TS",
                "Sci peak TS"
                ]
            ytitle = f"FERS C/S {gain}"
            for hist_name, xtitle in zip(hist_names, xtitles):
                output_name = hist_name.replace("hist_","")
                hist = infile.Get(hist_name)
                extraToDraw = ROOT.TPaveText(0.20, 0.85, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(f"correlation = {hist.GetCorrelationFactor():.3f}")
                DrawHistos([hist], "", TSmin if "Sci" in xtitle else TSCermin, TSmax if "Sci" in xtitle else TSCermax, xtitle, 0.0,2.0, ytitle,
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
                       outdir=outdir_plots, addOverflow=False, run_number=run_number_text,extraToDraw=[extraToDraw])
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/FERSvsDRS//CSratio_VS_DRSTS.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html


def main():
    output_html_DRSPeakTSVSEnergy = makeDRSPeakTSVSEnergyPlots()
    print(f"DRS Peak TS VS energy plots saved to {output_html_DRSPeakTSVSEnergy}")
    output_html_DRSPeakTSVSCSratio = makeDRSPeakTSVSCSratioPlots()
    print(f"DRS Peak TS VS C/S plots saved to {output_html_DRSPeakTSVSCSratio}")

if __name__ == "__main__":
    main()
