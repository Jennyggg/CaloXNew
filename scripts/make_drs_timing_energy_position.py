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
import matplotlib.pyplot as plt
import numpy as np

run_number = 1513
run_numbers_text = ["1501","1507","1511","1513","1513_2"]
positions = np.array([-168,-218,-268,-54.5,-400.3])
benergy = 40
run_number_text = "1501,1507,1511,1513"
plotdir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/plots/RunElectron90norm"
rootdir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/root/Run{}"
htmldir = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/html/RunElectron90norm"
GainCalibs = [("HG", False), ("LG", False), ("Mix", True)]
HE = (benergy >= 50)
TSmin = -80
TSmax = -20
TSCermin = -80
TSCermax = -50
def makeDRSPeakTSVSEnergyNormPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeakTS_VS_Energy"
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
                output_hist = None
                for run_text in run_numbers_text:
                    infile_name = rootdir.format(run_text)+"/drspeak_vs_fers.root"
                    print("open ",infile_name)
                    infile = ROOT.TFile(infile_name, "READ")
                    hist = infile.Get(hist_name)
                    print("Get ",hist_name)
                    if output_hist is None:
                        output_hist = hist.Clone()
                        output_hist.SetDirectory(0)
                        output_hist.Scale(1/output_hist.Integral()/len(run_numbers_text))
                        print("output_hist",output_hist)
                    else:
                        hist.Scale(1/hist.Integral()/len(run_numbers_text))
                        output_hist.Add(hist)
                        output_hist.SetDirectory(0)
                        print("output_hist",output_hist)
                    print("output_hist before closing",output_hist)
                    infile.Close()
                    print("output_hist after closing",output_hist)
                extraToDraw = ROOT.TPaveText(0.20, 0.85, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(f"correlation = {output_hist.GetCorrelationFactor():.3f}")
                DrawHistos([output_hist], "", TSmin if "Sci" in xtitle else TSCermin, TSmax if "Sci" in xtitle else TSCermax, xtitle, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"], ytitle,
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1e-4, zmax=output_hist.GetMaximum()*1.1, dologz=True,
                       outdir=outdir_plots, addOverflow=False, run_number=run_number_text,extraToDraw=[extraToDraw])
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/FERSvsDRS//EnergySum_VS_DRSTS.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html


def makeDRSPeakTSVSPositionPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeakTS_VS_Energy"
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
            DRS_names = [
                "DRSPeakTS_Cer_Quartz",
                "DRSPeakTS_Cer_Plastic",
                "DRSPeakTS_Sci"
            ]
            FERS_name = f"FERS_{cat}_{gain}"
            for hist_name, xtitle,DRS_name in zip(hist_names, xtitles,DRS_names):
                mean_X = []
                mean_X_err = []
                mean_Y = []
                mean_Y_err = []
                output_name = hist_name.replace("hist_","")
                for run_text in run_numbers_text:
                    infile_name = rootdir.format(run_text)+"/drspeak_vs_fers.root"
                    infile = ROOT.TFile(infile_name, "READ")
                    hist = infile.Get(hist_name)
                    hprojx = hist.ProjectionX("hprojx")
                    xmin = hprojx.GetMean() - 1.5*hprojx.GetRMS()
                    xmax = hprojx.GetMean() + 1.5*hprojx.GetRMS()
                    gausx = ROOT.TF1("gaus", "gaus", xmin, xmax)
                    hprojx.Fit(gausx, "RQ")
                    mean_X.append(gausx.GetParameter(1))
                    mean_X_err.append(gausx.GetParError(1))
                    hprojy = hist.ProjectionY("hprojy")
                    ymin = hprojy.GetMean() - 2*hprojy.GetRMS()
                    ymax = hprojy.GetMean() + 2*hprojy.GetRMS()
                    gausy = ROOT.TF1("gaus", "gaus", ymin, ymax)
                    hprojy.Fit(gausy, "RQ")
                    mean_Y.append(gausy.GetParameter(1))
                    mean_Y_err.append(gausy.GetParError(1))
                    infile.Close()
                mean_X = np.array(mean_X)
                mean_X_err = np.array(mean_X_err)
                mean_Y = np.array(mean_Y)
                mean_Y_err = np.array(mean_Y_err)
                if cat == "sci":
                    coeffs_x, cov_x = np.polyfit(
                        positions,
                        mean_X,
                        deg=1,
                        w=1.0 / mean_X_err,
                        cov=True
                        )
                    slope_x, intercept_x = coeffs_x
                    slope_err_x, intercept_err_x = np.sqrt(np.diag(cov_x))
                    pos_fit = np.linspace(positions.min(), positions.max(), 200)
                    x_fit = slope_x * pos_fit + intercept_x
                    plt.figure(figsize=(7, 5))
                    plt.errorbar(
                        positions,
                        mean_X,
                        yerr=mean_X_err,
                        fmt='o',
                        capsize=3,
                        label="90 degree e"
                    )
                    plt.plot(
                        pos_fit,
                        x_fit,
                        label=f"Fit: TSpeak = ({slope_x:.3f}±{slope_err_x:.3f})x + ({intercept_x:.3f}±{intercept_err_x:.3f})"
                        )
                    plt.xlabel("Table position (mm)")
                    plt.ylabel(xtitle)
                    plt.legend()
                    plt.grid(True)
                    plt.tight_layout()
                    plt.savefig(f"{outdir_plots}/{DRS_name}_vs_position.png", dpi=150)


                    coeffs_x, cov_x = np.polyfit(
                        positions,
                        mean_X/5,
                        deg=1,
                        w=1.0 / mean_X_err * 5,
                        cov=True
                        )
                    slope_x, intercept_x = coeffs_x
                    slope_err_x, intercept_err_x = np.sqrt(np.diag(cov_x))
                    x_fit = slope_x * pos_fit + intercept_x

                    plt.figure(figsize=(7, 5))
                    plt.errorbar(
                        positions,
                        mean_X/5,
                        yerr=mean_X_err/5,
                        fmt='o',
                        capsize=3,
                        label="90 degree e"
                    )
                    plt.plot(
                        pos_fit,
                        x_fit,
                        label=f"Fit: peak time = ({slope_x:.3f}±{slope_err_x:.3f})x + ({intercept_x:.3f}±{intercept_err_x:.3f})"
                        )
                    plt.xlabel("Table position (mm)")
                    plt.ylabel(xtitle.replace("TS","time")+" [ns]")
                    plt.legend()
                    plt.grid(True)
                    plt.tight_layout()
                    DRS_time_name = DRS_name.replace("TS","time")
                    plt.savefig(f"{outdir_plots}/{DRS_time_name}_vs_position.png", dpi=150)
                if DRS_name == "DRSPeakTS_Cer_Plastic":
                    coeffs_y, cov_y = np.polyfit(
                        positions,
                        mean_Y,
                        deg=1,
                        w=1.0 / mean_Y_err,
                        cov=True
                        )
                    slope_y, intercept_y = coeffs_y
                    slope_err_y, intercept_err_y = np.sqrt(np.diag(cov_y))
                    pos_fit = np.linspace(positions.min(), positions.max(), 200)
                    y_fit = slope_y * pos_fit + intercept_y
                    plt.figure(figsize=(7, 5))
                    plt.errorbar(
                        positions,
                        mean_Y,
                        yerr=mean_Y_err,
                        fmt='o',
                        capsize=3,
                        label="90 degree e"
                    )
                    plt.plot(
                        pos_fit,
                        y_fit,
                        label=f"Fit: Energy = ({slope_y:.3f}±{slope_err_y:.3f})x + ({intercept_y:.3f}±{intercept_err_y:.3f})"
                        )
                    plt.xlabel("Table position (mm)")
                    plt.ylabel(ytitle)
                    plt.legend()
                    plt.grid(True)
                    plt.tight_layout()
                    plt.savefig(f"{outdir_plots}/{FERS_name}_vs_position.png", dpi=150)


        
def main():
    output_html_DRSPeakTSVSEnergy = makeDRSPeakTSVSEnergyNormPlots()
    print(f"DRS Peak TS VS energy plots saved to {output_html_DRSPeakTSVSEnergy}")
    makeDRSPeakTSVSPositionPlots()

if __name__ == "__main__":
    main()
