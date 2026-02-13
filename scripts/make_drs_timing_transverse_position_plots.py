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
import math
auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
run_number = args.run

analysis = (CaloXAnalysisManager(args)
            .prepare()                   # Baseline and vectorization
            .apply_hole_veto(flag_only=True)
            )

fersboards = analysis.fersboards
DRSBoards = analysis.drsboards

benergy = analysis.beam_energy
run_number = analysis.run_number
paths = analysis.paths
rootdir = paths["root"]
plotdir = paths["plots"]
htmldir = paths["html"]

file_drschannels_bad = "data/drs/badchannels.json"
with open(file_drschannels_bad, "r") as f:
    drschannels_bad = json.load(f)

# rdf = analysis.get_rdf()
rdf = analysis.get_particle_analysis("muon")

TSmin = -80
TSmax = -20
TSCermin = -80
TSCermax = -60

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

def GetBinIntex(x,start,end,nbins):
    return math.ceil( (x-start)/(end-start) * nbins )

def checkDRSPeakTSvsPosition(rdf):


    h2_DRSPeakTS_x_quartz = ROOT.TH2D("h2_DRSPeakTS_x_quartz", "DRS Cer Quartz Peak TS VS x position;Cer Quartz Peak TS;x", 
                TSCermax - TSCermin, TSCermin, TSCermax, 
                4, -2, 2)
    h2_DRSPeakTS_x_plastic = ROOT.TH2D("h2_DRSPeakTS_x_plastic", "DRS Cer Plastic Peak TS VS x position;Cer Plastic Peak TS;x", 
                TSCermax - TSCermin, TSCermin, TSCermax, 
                4, -2, 2)
    h2_DRSPeakTS_x_sci = ROOT.TH2D("h2_DRSPeakTS_x_sci", "DRS Sci Peak TS VS x position;Sci Peak TS;x", 
                TSmax - TSmin, TSmin, TSmax, 
                4, -2, 2)
    h2_DRSPeakTS_y_quartz = ROOT.TH2D("h2_DRSPeakTS_y_quartz", "DRS Cer Quartz Peak TS VS y position;Cer Quartz Peak TS;y", 
                TSCermax - TSCermin, TSCermin, TSCermax, 
                16, -2, 2)
    h2_DRSPeakTS_y_plastic = ROOT.TH2D("h2_DRSPeakTS_y_plastic", "DRS Cer Plastic Peak TS VS y position;Cer Plastic Peak TS;y", 
                TSCermax - TSCermin, TSCermin, TSCermax, 
                16, -2, 2)
    h2_DRSPeakTS_y_sci = ROOT.TH2D("h2_DRSPeakTS_y_sci", "DRS Sci Peak TS VS y position;Sci Peak TS;y", 
                TSmax - TSmin, TSmin, TSmax, 
                16, -2, 2)
    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        # if board_no > 3:
        #    continue
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)
            if i_tower_x < -2 or i_tower_x > 2: continue
            if i_tower_y < -2 or i_tower_y > 2: continue
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
                    TSmaxtmp = TSCermax
                    TSmintmp = TSCermin
                    thre_Sum = 1000
                else:
                    TSmaxtmp = TSmax
                    TSmintmp = TSmin
                    thre_Sum = 5000
                rdf = rdf.Define(f"{channelName}_{varsuffix}_hasSignal",f"(int)({channelName}_Sum > {thre_Sum}) * (int)({channelName}_{varsuffix} < {TSmaxtmp}) * (int)({channelName}_{varsuffix} > {TSmintmp}) * {channelName}_{varsuffix}")
                h1_DRSPeakTS = rdf.Histo1D((
                    f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Peak TS for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var};Peak TS;Counts",
                    TSmaxtmp - TSmintmp, TSmintmp, TSmaxtmp),
                    f"{channelName}_{varsuffix}_hasSignal"
                )

                if var == "Cer":
                    if chan_DRS.isQuartz:
                        for ibin in range(1, TSmaxtmp - TSmintmp + 1):
                            h2_DRSPeakTS_x_quartz.SetBinContent(ibin,GetBinIntex(i_tower_x,-2,2,4),h2_DRSPeakTS_x_quartz.GetBinContent(ibin,GetBinIntex(i_tower_x,-2,2,4))+h1_DRSPeakTS.GetBinContent(ibin))
                            h2_DRSPeakTS_x_quartz.SetBinError(ibin,GetBinIntex(i_tower_x,-2,2,4),math.sqrt(h2_DRSPeakTS_x_quartz.GetBinError(ibin,GetBinIntex(i_tower_x,-2,2,4))**2+h1_DRSPeakTS.GetBinError(ibin)**2))
                            h2_DRSPeakTS_y_quartz.SetBinContent(ibin,GetBinIntex(i_tower_y,-2,2,16),h2_DRSPeakTS_y_quartz.GetBinContent(ibin,GetBinIntex(i_tower_y,-2,2,16))+h1_DRSPeakTS.GetBinContent(ibin))
                            h2_DRSPeakTS_y_quartz.SetBinError(ibin,GetBinIntex(i_tower_y,-2,2,16),math.sqrt(h2_DRSPeakTS_y_quartz.GetBinError(ibin,GetBinIntex(i_tower_y,-2,2,16))**2+h1_DRSPeakTS.GetBinError(ibin)**2))
                    else:
                        for ibin in range(1, TSmaxtmp - TSmintmp + 1):
                            h2_DRSPeakTS_x_plastic.SetBinContent(ibin,GetBinIntex(i_tower_x,-2,2,4),h2_DRSPeakTS_x_plastic.GetBinContent(ibin,GetBinIntex(i_tower_x,-2,2,4))+h1_DRSPeakTS.GetBinContent(ibin))
                            h2_DRSPeakTS_x_plastic.SetBinError(ibin,GetBinIntex(i_tower_x,-2,2,4),math.sqrt(h2_DRSPeakTS_x_plastic.GetBinError(ibin,GetBinIntex(i_tower_x,-2,2,4))**2+h1_DRSPeakTS.GetBinError(ibin)**2))
                            h2_DRSPeakTS_y_plastic.SetBinContent(ibin,GetBinIntex(i_tower_y,-2,2,16),h2_DRSPeakTS_y_plastic.GetBinContent(ibin,GetBinIntex(i_tower_y,-2,2,16))+h1_DRSPeakTS.GetBinContent(ibin))
                            h2_DRSPeakTS_y_plastic.SetBinError(ibin,GetBinIntex(i_tower_y,-2,2,16),math.sqrt(h2_DRSPeakTS_y_plastic.GetBinError(ibin,GetBinIntex(i_tower_y,-2,2,16))**2+h1_DRSPeakTS.GetBinError(ibin)**2))
                else:
                    for ibin in range(1, TSmaxtmp - TSmintmp + 1):
                        h2_DRSPeakTS_x_sci.SetBinContent(ibin,GetBinIntex(i_tower_x,-2,2,4),h2_DRSPeakTS_x_sci.GetBinContent(ibin,GetBinIntex(i_tower_x,-2,2,4))+h1_DRSPeakTS.GetBinContent(ibin))
                        h2_DRSPeakTS_x_sci.SetBinError(ibin,GetBinIntex(i_tower_x,-2,2,4),math.sqrt(h2_DRSPeakTS_x_sci.GetBinError(ibin,GetBinIntex(i_tower_x,-2,2,4))**2+h1_DRSPeakTS.GetBinError(ibin)**2))
                        h2_DRSPeakTS_y_sci.SetBinContent(ibin,GetBinIntex(i_tower_y,-2,2,16),h2_DRSPeakTS_y_sci.GetBinContent(ibin,GetBinIntex(i_tower_y,-2,2,16))+h1_DRSPeakTS.GetBinContent(ibin))
                        h2_DRSPeakTS_y_sci.SetBinError(ibin,GetBinIntex(i_tower_y,-2,2,16),math.sqrt(h2_DRSPeakTS_y_sci.GetBinError(ibin,GetBinIntex(i_tower_y,-2,2,16))**2+h1_DRSPeakTS.GetBinError(ibin)**2))
            if len(channelNames) < 2:
                print(
                    f"Warning: Not enough good channels found for Board{board_no}, Tower({sTowerX}, {sTowerY})")
                continue


    return h2_DRSPeakTS_x_quartz,h2_DRSPeakTS_x_plastic, h2_DRSPeakTS_x_sci, h2_DRSPeakTS_y_quartz,h2_DRSPeakTS_y_plastic,h2_DRSPeakTS_y_sci


def makeDRSPeakTSvsPositionPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRS_vs_position"
    infile_name = f"{rootdir}/drspeakts_vs_position.root"
    infile = ROOT.TFile(infile_name, "READ")
    for hist2d_name in ["h2_DRSPeakTS_x_quartz","h2_DRSPeakTS_x_plastic","h2_DRSPeakTS_x_sci","h2_DRSPeakTS_y_quartz","h2_DRSPeakTS_y_plastic","h2_DRSPeakTS_y_sci"]:
        hist2d = infile.Get(hist2d_name)
        output_name = hist2d_name
        plots.append(output_name + ".png")
        if not hist2d:
            print(
                f"Warning: Histogram {hist2d_name} not found in {infile_name}")
            continue
        xtitle = hist2d.GetXaxis().GetTitle()
        ytitle = hist2d.GetYaxis().GetTitle()
        DrawHistos([hist2d], "", TSmin if "Sci" in xtitle else TSCermin, TSmax if "Sci" in xtitle else TSCermax, xtitle, -2, 2, ytitle,
                    output_name,
                    dology=False, addOverflow=False, addUnderflow=False, doth2=True, drawoptions="COLZ", zmin=1, zmax=1e4, dologz=True,
                    outdir=outdir_plots, run_number=run_number)

        plots.append(output_name + ".png")
        print(outdir_plots+ "/" +output_name + ".png saved")
    output_html = f"{htmldir}/DRSPosition/DRS_VS_position.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html




def main():
    makeHists = True
    makePlots = True

    if makeHists:
        global rdf

        rdf = calibrateDRSPeakTS(rdf, run_number, DRSBoards,
                                 TSminDRS=450, TSmaxDRS=550, threshold=200.0)

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

        h2_DRSPeakTS_x_quartz,h2_DRSPeakTS_x_plastic, h2_DRSPeakTS_x_sci, h2_DRSPeakTS_y_quartz,h2_DRSPeakTS_y_plastic,h2_DRSPeakTS_y_sci = checkDRSPeakTSvsPosition(rdf)
        outfile_DRSPeakTSVSPosition = ROOT.TFile(
            f"{rootdir}/drspeakts_vs_position.root", "RECREATE")
        for hist in [h2_DRSPeakTS_x_quartz,h2_DRSPeakTS_x_plastic, h2_DRSPeakTS_x_sci, h2_DRSPeakTS_y_quartz,h2_DRSPeakTS_y_plastic,h2_DRSPeakTS_y_sci]:
            hist.SetDirectory(outfile_DRSPeakTSVSPosition)
            hist.Write()
        outfile_DRSPeakTSVSPosition.Close()
        print(f"{rootdir}/drspeakts_vs_position.root created")
    if makePlots:
        output_html_DRSPeakTSVSPosition = makeDRSPeakTSvsPositionPlots()

if __name__ == "__main__":
    main()