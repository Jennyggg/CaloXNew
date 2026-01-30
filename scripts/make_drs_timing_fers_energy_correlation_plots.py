import ROOT
import json
from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from configs.plot_config import getRangesForFERSEnergySums
from variables.drs import calibrateDRSPeakTS
from utils.html_generator import generate_html
from utils.parser import get_args
from plotting.my_function import DrawHistos, LHistos2Hist
from utils.timing import auto_timer
from utils.root_setup import setup_root
from utils.utils import number_to_string

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
run_number = args.run

analysis = (CaloXAnalysisManager(args)
            .prepare()                   # Baseline and vectorization
            .calibrate_fers()
            .apply_hole_veto(flag_only=True)
            )


GainCalibs = [("HG", False), ("LG", False), ("Mix", True)]

# calculate energy sums
for gain, calib in GainCalibs:
    analysis = analysis.define_physics_variables(
        gain=gain, calib=calib, pdsub=True)

fersboards = analysis.fersboards
DRSBoards = analysis.drsboards

benergy = analysis.beam_energy
run_number = analysis.run_number
paths = analysis.paths
rootdir = paths["root"]
plotdir = paths["plots"]
htmldir = paths["html"]

doPerBoardPlots = False
HE = (benergy >= 50)

file_drschannels_bad = "data/drs/badchannels.json"
with open(file_drschannels_bad, "r") as f:
    drschannels_bad = json.load(f)

rdf = analysis.get_particle_analysis("pion")

#TSmin = -90
#TSmax = -10
#TSCermin = -70
#TSCermax = -50



TSmin = -200
TSmax = 10
TSCermin = -400
TSCermax = 0

varsuffix = "DiffRelPeakTS_US"

def GetMean(hist, useMode=True):
    if useMode:
        # peak position of peakTS
        bin_max = hist.GetMaximumBin()
        mean = hist.GetBinCenter(bin_max)
    else:
        # average position of peakTS
        mean = hist.GetMean()
    return mean

def checkDRSPeakTSVSEnergy(rdf):

    hists_DRSPeakTS_vs_FERS_EnergySum = []
    channelnames_quartz = []
    channelnames_plastic = []
    channelnames_sci = []

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        # if board_no > 3:
        #    continue
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            channelNames = {}
            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_DRS is None:
                    print(
                        f"Warning: DRS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                channelName = chan_DRS.get_channel_name(blsub=False)
                if channelName in drschannels_bad:
                    print(
                        f"Warning: DRS Channel {channelName} is marked as bad channel, skipping...")
                    continue
                channelNames[var] = channelName

                if var == "Cer":
                    if chan_DRS.isQuartz:
                        channelnames_quartz.append(
                            f"{channelName}_{varsuffix}")
                    else:
                        channelnames_plastic.append(
                            f"{channelName}_{varsuffix}")
                else:
                    channelnames_sci.append(f"{channelName}_{varsuffix}")

            if len(channelNames) < 2:
                print(
                    f"Warning: Not enough good channels found for Board{board_no}, Tower({sTowerX}, {sTowerY})")
                continue


    # average of quartz, plastic and sci channels
    rdf = rdf.Define("Cer_Quartz_AvgPeakTS",
                     f"({'+'.join(channelnames_quartz)})/{len(channelnames_quartz)}")
    rdf = rdf.Define("Cer_Plastic_AvgPeakTS",
                     f"({'+'.join(channelnames_plastic)})/{len(channelnames_plastic)}")
    rdf = rdf.Define("Sci_AvgPeakTS",
                     f"({'+'.join(channelnames_sci)})/{len(channelnames_sci)}")

    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)

        for cat in ["cer", "sci"]:
            # per-event sum
            varname = fersboards.get_energy_sum_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)

            h2_DRSPeak_Cer_Quartz_VS_FERS = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Cer Quartz;Cer Quartz Peak TS;{cat} {gain} Energy",
                int((TSCermax - TSCermin)/5), TSCermin, TSCermax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Cer_Quartz_AvgPeakTS",
                varname
                )

            h2_DRSPeak_Cer_Plastic_VS_FERS = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Cer Plastic;Cer Plastic Peak TS;{cat} {gain} Energy",
                int((TSCermax - TSCermin)/5), TSCermin, TSmax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Cer_Plastic_AvgPeakTS",
                varname
                )

            h2_DRSPeak_Sci_VS_FERS = rdf.Histo2D((
                f"hist_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Sci ;Sci Peak TS;{cat} {gain} Energy",
                int((TSmax - TSmin)/5), TSmin, TSmax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Sci_AvgPeakTS",
                varname
                )

            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Cer_Quartz_VS_FERS)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Cer_Plastic_VS_FERS)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Sci_VS_FERS)


    return hists_DRSPeakTS_vs_FERS_EnergySum

def checkDRSPeakTSVSCSratio(rdf):

    hists_DRSPeakTS_vs_FERS_CSratio = []
    channelnames_quartz = []
    channelnames_plastic = []
    channelnames_sci = []

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        # if board_no > 3:
        #    continue
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            channelNames = {}
            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_DRS is None:
                    print(
                        f"Warning: DRS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                channelName = chan_DRS.get_channel_name(blsub=False)
                if channelName in drschannels_bad:
                    print(
                        f"Warning: DRS Channel {channelName} is marked as bad channel, skipping...")
                    continue
                channelNames[var] = channelName

                if var == "Cer":
                    if chan_DRS.isQuartz:
                        channelnames_quartz.append(
                            f"{channelName}_{varsuffix}")
                    else:
                        channelnames_plastic.append(
                            f"{channelName}_{varsuffix}")
                else:
                    channelnames_sci.append(f"{channelName}_{varsuffix}")

            if len(channelNames) < 2:
                print(
                    f"Warning: Not enough good channels found for Board{board_no}, Tower({sTowerX}, {sTowerY})")
                continue


    # average of quartz, plastic and sci channels
    rdf = rdf.Define("Cer_Quartz_AvgPeakTS",
                     f"({'+'.join(channelnames_quartz)})/{len(channelnames_quartz)}")
    rdf = rdf.Define("Cer_Plastic_AvgPeakTS",
                     f"({'+'.join(channelnames_plastic)})/{len(channelnames_plastic)}")
    rdf = rdf.Define("Sci_AvgPeakTS",
                     f"({'+'.join(channelnames_sci)})/{len(channelnames_sci)}")

    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)
        varname_cer = fersboards.get_energy_sum_name(
                gain=gain, isCer=True, pdsub=True, calib=calib)
        varname_sci = fersboards.get_energy_sum_name(
                gain=gain, isCer=False, pdsub=True, calib=calib)
        varname_ratio = f"{varname_cer}_{varname_sci}_ratio"
        rdf = rdf.Define(varname_ratio,
                     f"{varname_cer}/{varname_sci}")
        h2_DRSPeak_Cer_Quartz_VS_FERSratio = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_Quartz_VS_FERS_CSratio_{gain}",
                f"FERS C/S {gain} VS DRS Peak TS - Cer Quartz;Cer Quartz Peak TS;C/S {gain}",
                int((TSCermax - TSCermin)/5), TSCermin, TSCermax,
                100, 0.0, 1.0),
                "Cer_Quartz_AvgPeakTS",
                varname_ratio
                )

        h2_DRSPeak_Cer_Plastic_VS_FERSratio = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_Plastic_VS_FERS_CSratio_{gain}",
                f"FERS C/S {gain} VS DRS Peak TS - Cer Plastic;Cer Plastic Peak TS;C/S {gain}",
                int((TSCermax - TSCermin)/5), TSCermin, TSmax,
                100, 0.0,1.0),
                "Cer_Plastic_AvgPeakTS",
                varname_ratio
                )

        h2_DRSPeak_Sci_VS_FERSratio = rdf.Histo2D((
                f"hist_DRSPeakTS_Sci_VS_FERS_CSratio_{gain}",
                f"FERS C/S {gain} VS DRS Peak TS - Sci ;Sci Peak TS;C/S {gain} Energy",
                int((TSmax - TSmin)/5), TSmin, TSmax,
                100, 0.0,1.0),
                "Sci_AvgPeakTS",
                varname_ratio
                )

        hists_DRSPeakTS_vs_FERS_CSratio.append(h2_DRSPeak_Cer_Quartz_VS_FERSratio)
        hists_DRSPeakTS_vs_FERS_CSratio.append(h2_DRSPeak_Cer_Plastic_VS_FERSratio)
        hists_DRSPeakTS_vs_FERS_CSratio.append(h2_DRSPeak_Sci_VS_FERSratio)


    return hists_DRSPeakTS_vs_FERS_CSratio


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
                       outdir=outdir_plots, addOverflow=False, run_number=run_number,extraToDraw=[extraToDraw])
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
                DrawHistos([hist], "", TSmin if "Sci" in xtitle else TSCermin, TSmax if "Sci" in xtitle else TSCermax, xtitle, 0.0,1.0, ytitle,
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
                       outdir=outdir_plots, addOverflow=False, run_number=run_number,extraToDraw=[extraToDraw])
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/FERSvsDRS//CSratio_VS_DRSTS.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html


def main():
    makeHists = True
    makePlots = True

    if makeHists:
        global rdf

        rdf = calibrateDRSPeakTS(rdf, run_number, DRSBoards,
                                 TSminDRS=0, TSmaxDRS=1000, threshold=9.0)

        rdf_prefilterMCP1 = rdf
        map_mcp_channels = get_mcp_channels(run_number)

        condition = f"{map_mcp_channels['US'][0]}_RelPeakTS > -350 && {map_mcp_channels['US'][0]}_RelPeakTS < -100"
        condition += f" && {map_mcp_channels['US'][0]}_PeakTS > 500 && {map_mcp_channels['US'][0]}_PeakTS < 600"
        condition += f" && {map_mcp_channels['US'][0]}_Peak < -300.0"
        rdf_prefilterMCP2 = rdf_prefilterMCP1.Filter(condition,
                                                     "Pre-filter on MCP US channel 0 Peak TS")

        rdf_prefilterMCP2 = rdf_prefilterMCP1.Define(
            "MCP0_DeltaRelPeakTS", f"{map_mcp_channels['DS'][0]}_RelPeakTS - {map_mcp_channels['US'][0]}_RelPeakTS")
        rdf = rdf_prefilterMCP2

        hists_DRSTS_Energy = checkDRSPeakTSVSEnergy(rdf)
        hists_DRSTS_CSratio = checkDRSPeakTSVSCSratio(rdf)
        outfile_DRSPeakTS_FERS = ROOT.TFile(
            f"{rootdir}/drspeak_vs_fers.root", "RECREATE")
        for hist in hists_DRSTS_Energy:
            hist.SetDirectory(outfile_DRSPeakTS_FERS)
            hist.Write()
        outfile_DRSPeakTS_CSratio = ROOT.TFile(
            f"{rootdir}/drspeak_vs_fers_csratio.root", "RECREATE")
        for hist in hists_DRSTS_CSratio:
            hist.SetDirectory(outfile_DRSPeakTS_CSratio)
            hist.Write()
    if makePlots:
        output_html_DRSPeakTSVSEnergy = makeDRSPeakTSVSEnergyPlots()
        print(f"DRS Peak TS VS energy plots saved to {output_html_DRSPeakTSVSEnergy}")
        output_html_DRSPeakTSVSCSratio = makeDRSPeakTSVSCSratioPlots()
        print(f"DRS Peak TS VS C/S plots saved to {output_html_DRSPeakTSVSCSratio}")

if __name__ == "__main__":
    main()

