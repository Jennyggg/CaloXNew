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
from scipy.optimize import curve_fit
import numpy as np
import json
import os

benergies = [30,40,60,100,120,160]
labels = ["30 GeV","40 GeV","60 GeV","100 GeV","120 GeV","160 GeV"]
run_numbers = [1439,1438,1437,1431,1433,1434]
colors = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kMagenta,ROOT.kOrange,ROOT.kCyan]
GainCalibs = [("HG", False), ("LG", False), ("Mix", True)]
plotdir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/plots/RunPionsCenter"
htmldir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/html/RunPionsCenter"
jsondir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/root/RunPionsCenter"

TSmin = -80
TSmax = -20
TSCermin = -80
TSCermax = -50
TSAnchor = {
    "Cer_Quartz": -63.06,
    "Cer_Plastic": -60.22,
    "Sci": -45.05
}

def fitDRSPeakTSVSEnergyProfile():
    plots = []
    hist_dir = {}
    fitpoints_dir = {}
    def func(x, m, c):
        return m * x + c
    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            hist_dir[f"prof_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}"] = []
            hist_dir[f"prof_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}"] = []
            hist_dir[f"prof_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}"] = []
            fitpoints_dir[f"fit_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}"] = {
                "DRSPeakTS": [],
                "response": [],
                "error": []
            }
            fitpoints_dir[f"fit_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}"] = {
                "DRSPeakTS": [],
                "response": [],
                "error": []
            }
            fitpoints_dir[f"fit_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}"] = {
                "DRSPeakTS": [],
                "response": [],
                "error": []
            }
    outdir_plots = f"{plotdir}/DRSPeakTS_VS_EnergyProfiles"
    for benergy,run_number in zip(benergies,run_numbers):
        rootdir = f"/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/root/Run{run_number}"
        infile_name = f"{rootdir}/drspeak_vs_fers.root"
        infile = ROOT.TFile(infile_name, "READ")
        for key in hist_dir.keys():
            prof = infile.Get(key)
            prof.Scale(1.0/benergy)
            hist2d = infile.Get(key.replace("prof_","hist_"))
            proj = hist2d.ProjectionX()
            mean = proj.GetMean()
            std = proj.GetRMS()
            ibin_min = prof.GetXaxis().FindBin(mean-1.5*std)
            ibin_max = prof.GetXaxis().FindBin(mean+1.5*std)

            for ibin in range(ibin_min,ibin_max+1):
                fitpoints_dir[key.replace("prof_","fit_")]["DRSPeakTS"].append(prof.GetBinCenter(ibin))
                fitpoints_dir[key.replace("prof_","fit_")]["response"].append(prof.GetBinContent(ibin))
                fitpoints_dir[key.replace("prof_","fit_")]["error"].append(prof.GetBinError(ibin))

            prof.Fit("pol1","R","",mean-1.5*std,mean+1.5*std)
            prof.SetDirectory(0)
            hist_dir[key].append(prof)
    
    for key in fitpoints_dir.keys():
        TS = np.array(fitpoints_dir[key]["DRSPeakTS"])
        response = np.array(fitpoints_dir[key]["response"])
        error = np.array(fitpoints_dir[key]["error"])
        popt, pcov = curve_fit(func, TS, response, sigma=error, absolute_sigma=True)
        fitpoints_dir[key]["slope"] = popt[0]
        fitpoints_dir[key]["intercept"] = popt[1]
        fitpoints_dir[key]["cov"] = pcov

    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            prof_names = [
                f"prof_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}",
                f"prof_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}",
                f"prof_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}",
                ]
            xtitles = [
                "Cer Quartz peak TS",
                "Cer Plastic peak TS",
                "Sci peak TS"
                ]
            ytitle = f"FERS Response {cat} {gain}"

            for prof_name, xtitle in zip(prof_names, xtitles):
                output_name = prof_name.replace("prof_","")+"_mergefit"
                merge_fit = ROOT.TF1("merge_fit","[0] + [1]*x", TSmin if "Sci" in xtitle else TSCermin, TSmax if "Sci" in xtitle else TSCermax)
                merge_fit.SetParameters(fitpoints_dir[prof_name.replace("prof_","fit_")]["intercept"], fitpoints_dir[prof_name.replace("prof_","fit_")]["slope"])
                merge_fit.SetLineColor(ROOT.kRed)
                merge_fit.SetLineWidth(2)
                extraToDraw = ROOT.TPaveText(0.20, 0.15, 0.60, 0.2, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(f"response = {fitpoints_dir[prof_name.replace('prof_','fit_')]['slope']:.3e} TS + {fitpoints_dir[prof_name.replace('prof_','fit_')]['intercept']:.3e}")

                DrawHistos(hist_dir[prof_name], labels, TSmin if "Sci" in xtitle else TSCermin, TSmax if "Sci" in xtitle else TSCermax, xtitle, 0, 2, ytitle,
                       output_name,
                       dology=False, mycolors=colors, addOverflow=False, addUnderflow=False,
                       outdir=outdir_plots, run_number="Pion runs",extraToDraw=[extraToDraw,merge_fit])
            plots.append(output_name + ".png")
    output_html = f"{htmldir}/DRSPeakTS_VS_EnergyProfiles//EnergySum_VS_DRSTS_profile_merge.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return fitpoints_dir,output_html


def makeDRSPeakTSVSEnergyProfilePlots():
    plots = []
    hist_dir = {}
    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            hist_dir[f"prof_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}"] = []
            hist_dir[f"prof_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}"] = []
            hist_dir[f"prof_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}"] = []
    outdir_plots = f"{plotdir}/DRSPeakTS_VS_EnergyProfiles"
    for benergy,run_number in zip(benergies,run_numbers):
        rootdir = f"/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/root/Run{run_number}"
        infile_name = f"{rootdir}/drspeak_vs_fers.root"
        infile = ROOT.TFile(infile_name, "READ")
        for key in hist_dir.keys():
            prof = infile.Get(key)
            prof.Scale(1.0/benergy)
            hist2d = infile.Get(key.replace("prof_","hist_"))
            proj = hist2d.ProjectionX()
            mean = proj.GetMean()
            std = proj.GetRMS()
            prof.Fit("pol1","R","",mean-1.5*std,mean+1.5*std)
            prof.SetDirectory(0)
            hist_dir[key].append(prof)
    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            prof_names = [
                f"prof_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}",
                f"prof_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}",
                f"prof_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}",
                ]
            xtitles = [
                "Cer Quartz peak TS",
                "Cer Plastic peak TS",
                "Sci peak TS"
                ]
            ytitle = f"FERS Response {cat} {gain}"

            for prof_name, xtitle in zip(prof_names, xtitles):
                output_name = prof_name.replace("prof_","")
                print(prof_name,hist_dir[prof_name])
                

                DrawHistos(hist_dir[prof_name], labels, TSmin if "Sci" in xtitle else TSCermin, TSmax if "Sci" in xtitle else TSCermax, xtitle, 0, 2 if gain == "Mix" else None, ytitle,
                       output_name,
                       dology=False, mycolors=colors, addOverflow=False, addUnderflow=False,
                       outdir=outdir_plots, run_number="Pion runs")
            plots.append(output_name + ".png")
    output_html = f"{htmldir}/DRSPeakTS_VS_EnergyProfiles//EnergySum_VS_DRSTS_profile.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html




def main():
    output_html_DRSPeakTSVSEnergyProfile = makeDRSPeakTSVSEnergyProfilePlots()
    print(f"DRS Peak TS VS response profile plots saved to {output_html_DRSPeakTSVSEnergyProfile}")
    fit, output_html_DRSPeakTSVSEnergyProfileMerge = fitDRSPeakTSVSEnergyProfile()

    fit_write = {}
    for key in fit.keys():
        fit_write[key] = {
            "slope": fit[key]["slope"],
            "intercept": fit[key]["intercept"],
        }
    os.makedirs(jsondir, exist_ok=True)
    with open(f"{jsondir}/fit_response_TS.json", "w") as json_file:
        json.dump(fit_write, json_file, indent=4)


if __name__ == "__main__":
    main()
