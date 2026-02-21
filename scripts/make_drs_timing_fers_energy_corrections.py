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
chi_DR = 0.12
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
corr_json = "/lustre/work/jweijie/CaloXNew/CaloXDataAnalysis/results/root/RunPionsCenter/fit_response_TS.json"
with open(corr_json, 'r') as file:
    corr_dic = json.load(file)

TSAnchor = {
    "Cer_Quartz": -63.06,
    "Cer_Plastic": -60.22,
    "Sci": -45.05
}

doPerBoardPlots = False
HE = (benergy >= 50)

file_drschannels_bad = "data/drs/badchannels.json"
with open(file_drschannels_bad, "r") as f:
    drschannels_bad = json.load(f)

rdf = analysis.get_particle_analysis("pion")

TSmin = -80
TSmax = -20
TSCermin = -80
TSCermax = -50


#TSmin = -200
#TSmax = 10
#TSCermin = -400
#TSCermax = 0

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

def checkDRSPeakTSVSEnergyCorr(rdf):

    hists_DRSPeakTS_vs_FERS_EnergySum = []
    proj_DRSPeakTS_vs_FERS_EnergySum = []
    channelnames_quartz = []
    channelnames_plastic = []
    channelnames_sci = []
    channelnames_weight_ave_quartz = []
    channelnames_weight_ave_plastic = []
    channelnames_weight_ave_sci = []

    channelnames_sum_quartz = []
    channelnames_sum_plastic = []
    channelnames_sum_sci = []

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
                        channelnames_weight_ave_quartz.append(
                            f"(int)({channelName}_Sum > 1800) * (int)({channelName}_{varsuffix} < -60) * (int)({channelName}_{varsuffix} > -80) * {channelName}_{varsuffix} * {channelName}_Sum"
                        )
                        channelnames_sum_quartz.append(
                            f"(int)({channelName}_Sum > 1800) * (int)({channelName}_{varsuffix} < -60) * (int)({channelName}_{varsuffix} > -80) * {channelName}_Sum"
                        )
                    else:
                        channelnames_plastic.append(
                            f"{channelName}_{varsuffix}")
                        channelnames_weight_ave_plastic.append(
                            f"(int)({channelName}_Sum > 2500) * (int)({channelName}_{varsuffix} < -50) * (int)({channelName}_{varsuffix} > -80) * {channelName}_{varsuffix} * {channelName}_Sum"
                        )
                        channelnames_sum_plastic.append(
                            f"(int)({channelName}_Sum > 2500) * (int)({channelName}_{varsuffix} < -50) * (int)({channelName}_{varsuffix} > -80) * {channelName}_Sum"
                        )
                else:
                    channelnames_sci.append(f"{channelName}_{varsuffix}")
                    channelnames_weight_ave_sci.append(f"(int)({channelName}_Sum>8000) * (int)({channelName}_{varsuffix} < -20) * (int)({channelName}_{varsuffix} > -80) * {channelName}_{varsuffix} * {channelName}_Sum")
                    channelnames_sum_sci.append(
                            f"(int)({channelName}_Sum >8000 ) * (int)({channelName}_{varsuffix} < -20) * (int)({channelName}_{varsuffix} > -80) * {channelName}_Sum"
                        )

            if len(channelNames) < 2:
                print(
                    f"Warning: Not enough good channels found for Board{board_no}, Tower({sTowerX}, {sTowerY})")
                continue


    # average of quartz, plastic and sci channels
    rdf = rdf.Define("Cer_Quartz_AvgPeakTS",
                     f"({'+'.join(channelnames_weight_ave_quartz)})/({'+'.join(channelnames_sum_quartz)})")
    rdf = rdf.Define("Cer_Plastic_AvgPeakTS",
                     f"({'+'.join(channelnames_weight_ave_plastic)})/({'+'.join(channelnames_sum_plastic)})")
    rdf = rdf.Define("Sci_AvgPeakTS",
                     f"({'+'.join(channelnames_weight_ave_sci)})/({'+'.join(channelnames_sum_sci)})")

    
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)

        for cat in ["cer", "sci"]:
            # per-event sum
            varname = fersboards.get_energy_sum_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
            slope = corr_dic[f"fit_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}"]["slope"]
            intercept = corr_dic[f"fit_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}"]["intercept"]
            T0 = TSAnchor["Cer_Quartz"]
            print(f"fit_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}")
            print("slope",slope,"intercept",intercept,"T0",T0)
            print("define ",f"{varname}*{slope*T0+intercept}/({slope}*Cer_Quartz_AvgPeakTS{'' if intercept<0 else '+'}{intercept})")
            rdf = rdf.Define(f"{varname}_Cer_Quartz_corrfactor",f"{slope*T0+intercept}/({slope}*Cer_Quartz_AvgPeakTS{'' if intercept<0 else '+'}{intercept})")
            rdf = rdf.Define(f"{varname}_corr_Cer_Quartz",f"{varname}*{varname}_Cer_Quartz_corrfactor")
            h2_DRSPeak_Cer_Quartz_VS_FERS = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Cer Quartz;Cer Quartz Peak TS;{cat} {gain} Energy",
                TSCermax - TSCermin, TSCermin, TSCermax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Cer_Quartz_AvgPeakTS",
                varname
                )
            print()
            proj_DRSPeak_Cer_Quartz_VS_FERS = h2_DRSPeak_Cer_Quartz_VS_FERS.ProjectionY(f"proj_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}",1,TSCermax - TSCermin)
            h2_DRSPeak_Cer_Quartz_VS_FERS_corr = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}_corr",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Cer Quartz;Cer Quartz Peak TS;{cat} {gain} Energy Corr",
                TSCermax - TSCermin, TSCermin, TSCermax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Cer_Quartz_AvgPeakTS",
                f"{varname}_corr_Cer_Quartz"
                )
            proj_DRSPeak_Cer_Quartz_VS_FERS_corr = h2_DRSPeak_Cer_Quartz_VS_FERS_corr.ProjectionY(f"proj_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}_corr",1,TSCermax - TSCermin)
            proj_DRSPeak_Cer_Quartz_VS_FERS_corr.SetLineColor(ROOT.kRed)
            proj_DRSPeak_Cer_Quartz_VS_FERS_corr.SetMarkerColor(ROOT.kRed)
            
            
            slope = corr_dic[f"fit_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}"]["slope"]
            intercept = corr_dic[f"fit_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}"]["intercept"]
            T0 = TSAnchor["Cer_Plastic"]
            print(f"fit_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}")
            print("slope",slope,"intercept",intercept,"T0",T0)
            print("define ",f"{varname}*{slope*T0+intercept}/({slope}*Cer_Plastic_AvgPeakTS{'' if intercept<0 else '+'}{intercept})")
            rdf = rdf.Define(f"{varname}_Cer_Plastic_corrfactor",f"{slope*T0+intercept}/({slope}*Cer_Plastic_AvgPeakTS{'' if intercept<0 else '+'}{intercept})")
            rdf = rdf.Define(f"{varname}_corr_Cer_Plastic",f"{varname}*{varname}_Cer_Plastic_corrfactor")
            h2_DRSPeak_Cer_Plastic_VS_FERS = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Cer Plastic;Cer Plastic Peak TS;{cat} {gain} Energy",
                TSCermax - TSCermin, TSCermin, TSCermax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Cer_Plastic_AvgPeakTS",
                varname
                )
            proj_DRSPeak_Cer_Plastic_VS_FERS = h2_DRSPeak_Cer_Plastic_VS_FERS.ProjectionY(f"proj_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}",1,TSCermax - TSCermin)
            h2_DRSPeak_Cer_Plastic_VS_FERS_corr = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}_corr",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Cer Plastic;Cer Plastic Peak TS;{cat} {gain} Energy Corr",
                TSCermax - TSCermin, TSCermin, TSCermax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Cer_Plastic_AvgPeakTS",
                f"{varname}_corr_Cer_Plastic"
                )
            proj_DRSPeak_Cer_Plastic_VS_FERS_corr = h2_DRSPeak_Cer_Plastic_VS_FERS_corr.ProjectionY(f"proj_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}_corr",1,TSCermax - TSCermin)
            proj_DRSPeak_Cer_Plastic_VS_FERS_corr.SetLineColor(ROOT.kRed)
            proj_DRSPeak_Cer_Plastic_VS_FERS_corr.SetMarkerColor(ROOT.kRed)

            
            slope = corr_dic[f"fit_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}"]["slope"]
            intercept = corr_dic[f"fit_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}"]["intercept"]
            T0 = TSAnchor["Sci"]
            print(f"fit_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}")
            print("slope",slope,"intercept",intercept,"T0",T0)
            print("define ",f"{varname}*{slope*T0+intercept}/({slope}*Sci_AvgPeakTS{'' if intercept<0 else '+'}{intercept})")
            rdf = rdf.Define(f"{varname}_Sci_corrfactor",f"{slope*T0+intercept}/({slope}*Sci_AvgPeakTS{'' if intercept<0 else '+'}{intercept})")
            rdf = rdf.Define(f"{varname}_corr_Sci",f"{varname} * {varname}_Sci_corrfactor")
            h2_DRSPeak_Sci_VS_FERS = rdf.Histo2D((
                f"hist_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Sci;Sci Peak TS;{cat} {gain} Energy",
                TSmax - TSmin, TSmin, TSmax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Sci_AvgPeakTS",
                varname
                )
            proj_DRSPeak_Sci_VS_FERS = h2_DRSPeak_Sci_VS_FERS.ProjectionY(f"proj_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}",1,TSmax - TSmin)
            h2_DRSPeak_Sci_VS_FERS_corr = rdf.Histo2D((
                f"hist_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}_corr",
                f"FERS Energy {cat} {gain} VS DRS Peak TS - Sci;Sci Peak TS;{cat} {gain} Energy Corr",
                TSmax - TSmin, TSmin, TSmax,
                100, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                "Sci_AvgPeakTS",
                f"{varname}_corr_Sci"
                )
            proj_DRSPeak_Sci_VS_FERS_corr = h2_DRSPeak_Sci_VS_FERS_corr.ProjectionY(f"proj_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}_corr",1,TSmax - TSmin)
            proj_DRSPeak_Sci_VS_FERS_corr.SetLineColor(ROOT.kRed)
            proj_DRSPeak_Sci_VS_FERS_corr.SetMarkerColor(ROOT.kRed)


            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Cer_Quartz_VS_FERS_corr)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Cer_Plastic_VS_FERS_corr)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Sci_VS_FERS_corr)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Cer_Quartz_VS_FERS)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Cer_Plastic_VS_FERS)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_DRSPeak_Sci_VS_FERS)
            proj_DRSPeakTS_vs_FERS_EnergySum.append(proj_DRSPeak_Cer_Quartz_VS_FERS)
            proj_DRSPeakTS_vs_FERS_EnergySum.append(proj_DRSPeak_Cer_Plastic_VS_FERS)
            proj_DRSPeakTS_vs_FERS_EnergySum.append(proj_DRSPeak_Sci_VS_FERS)
            proj_DRSPeakTS_vs_FERS_EnergySum.append(proj_DRSPeak_Cer_Quartz_VS_FERS_corr)
            proj_DRSPeakTS_vs_FERS_EnergySum.append(proj_DRSPeak_Cer_Plastic_VS_FERS_corr)
            proj_DRSPeakTS_vs_FERS_EnergySum.append(proj_DRSPeak_Sci_VS_FERS_corr)

        varname_cer = fersboards.get_energy_sum_name(
                gain=gain, isCer=True, pdsub=True, calib=calib)
        varname_sci = fersboards.get_energy_sum_name(
                gain=gain, isCer=False, pdsub=True, calib=calib)

        h2_FERS_Cer_VS_Sci_filter_DRSPeak_Cer_Quartz = rdf.Filter(f"Cer_Quartz_AvgPeakTS>{TSCermin} && Cer_Quartz_AvgPeakTS<{TSCermax}").Histo2D((
                f"hist_FERS_{gain}_Cer_VS_Sci_filter_DRSPeak_Cer_Quartz",
                f"FERS Energy {cat} {gain} Cer VS Sci filtered by DRS Peak TS - Cer Quartz; {gain} Sci Energy;{gain} Cer Energy",
                100, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
                100, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
                varname_sci,
                varname_cer
                )
        h2_FERS_Cer_VS_Sci_DRSPeak_Cer_Quartz_corr = rdf.Filter(f"Cer_Quartz_AvgPeakTS>{TSCermin} && Cer_Quartz_AvgPeakTS<{TSCermax}").Histo2D((
                f"hist_FERS_{gain}_Cer_VS_Sci_DRSPeak_Cer_Quartz_corr",
                f"FERS Energy {cat} {gain} Cer VS Sci corrected by DRS Peak TS - Cer Quartz; {gain} Sci Energy;{gain} Cer Energy",
                100, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
                100, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
                f"{varname_sci}_corr_Cer_Quartz",
                f"{varname_cer}_corr_Cer_Quartz"
                )

        h2_FERS_Cer_VS_Sci_filter_DRSPeak_Cer_Plastic = rdf.Filter(f"Cer_Plastic_AvgPeakTS>{TSCermin} && Cer_Plastic_AvgPeakTS<{TSCermax}").Histo2D((
                f"hist_FERS_{gain}_Cer_VS_Sci_filter_DRSPeak_Cer_Plastic",
                f"FERS Energy {cat} {gain} Cer VS Sci filtered by DRS Peak TS - Cer Plastic; {gain} Sci Energy;{gain} Cer Energy",
                100, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
                100, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
                varname_sci,
                varname_cer
                )

        h2_FERS_Cer_VS_Sci_DRSPeak_Cer_Plastic_corr = rdf.Filter(f"Cer_Plastic_AvgPeakTS>{TSCermin} && Cer_Plastic_AvgPeakTS<{TSCermax}").Histo2D((
                f"hist_FERS_{gain}_Cer_VS_Sci_DRSPeak_Cer_Plastic_corr",
                f"FERS Energy {cat} {gain} Cer VS Sci corrected by DRS Peak TS - Cer Plastic; {gain} Sci Energy;{gain} Cer Energy",
                100, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
                100, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
                f"{varname_sci}_corr_Cer_Plastic",
                f"{varname_cer}_corr_Cer_Plastic"
                )

        h2_FERS_Cer_VS_Sci_filter_DRSPeak_Sci = rdf.Filter(f"Sci_AvgPeakTS>{TSmin} && Sci_AvgPeakTS<{TSmax}").Histo2D((
                f"hist_FERS_{gain}_Cer_VS_Sci_filter_DRSPeak_Sci",
                f"FERS Energy {cat} {gain} Cer VS Sci filtered by DRS Peak TS - Sci; {gain} Sci Energy;{gain} Cer Energy",
                100, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
                100, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
                varname_sci,
                varname_cer
                )
        h2_FERS_Cer_VS_Sci_DRSPeak_Sci_corr = rdf.Filter(f"Sci_AvgPeakTS>{TSmin} && Sci_AvgPeakTS<{TSmax}").Histo2D((
                f"hist_FERS_{gain}_Cer_VS_Sci_DRSPeak_Sci_corr",
                f"FERS Energy {cat} {gain} Cer VS Sci corrected by DRS Peak TS - Sci; {gain} Sci Energy;{gain} Cer Energy",
                100, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
                100, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
                f"{varname_sci}_corr_Sci",
                f"{varname_cer}_corr_Sci"
                )
        if gain == "Mix":
            rdf = rdf.Define("Energy_DR",f"({varname_sci}-{chi_DR}*{varname_cer})/(1-{chi_DR})")
            rdf = rdf.Define("Energy_DR_corr_Cer_Plastic",f"({varname_sci}_corr_Cer_Plastic-{chi_DR}*{varname_cer}_corr_Cer_Plastic)/(1-{chi_DR})")
            rdf = rdf.Define("Energy_DR_corr_Cer_Quartz",f"({varname_sci}_corr_Cer_Quartz-{chi_DR}*{varname_cer}_corr_Cer_Quartz)/(1-{chi_DR})")
            rdf = rdf.Define("Energy_DR_corr_Sci",f"({varname_sci}_corr_Sci-{chi_DR}*{varname_cer}_corr_Sci)/(1-{chi_DR})")
            hist_FERS_DR_filter_DRSPeak_Cer_Plastic = rdf.Filter(f"Cer_Plastic_AvgPeakTS>{TSCermin} && Cer_Plastic_AvgPeakTS<{TSCermax}").Histo1D((
                "hist_FERS_DR_filter_DRSPeak_Cer_Plastic",
                f"Combined energy with DRS Cer Plastic signals;Combined Energy;Counts",
                100, config["xmin_total"]["Mix_sci"], config["xmax_total"]["Mix_sci"]),
                "Energy_DR"
                )
            hist_FERS_DR_DRSPeak_Cer_Plastic_corr = rdf.Filter(f"Cer_Plastic_AvgPeakTS>{TSCermin} && Cer_Plastic_AvgPeakTS<{TSCermax}").Histo1D((
                "hist_FERS_DR_DRSPeak_Cer_Plastic_corr",
                f"Combined energy corrected by DRS Cer Plastic signals;Combined Energy;Counts",
                100, config["xmin_total"]["Mix_sci"], config["xmax_total"]["Mix_sci"]),
                "Energy_DR_corr_Cer_Plastic"
                )

            hist_FERS_DR_filter_DRSPeak_Cer_Quartz = rdf.Filter(f"Cer_Quartz_AvgPeakTS>{TSCermin} && Cer_Quartz_AvgPeakTS<{TSCermax}").Histo1D((
                "hist_FERS_DR_filter_DRSPeak_Cer_Quartz",
                f"Combined energy with DRS Cer Quartz signals;Combined Energy;Counts",
                100, config["xmin_total"]["Mix_sci"], config["xmax_total"]["Mix_sci"]),
                "Energy_DR"
                )
            hist_FERS_DR_DRSPeak_Cer_Quartz_corr = rdf.Filter(f"Cer_Quartz_AvgPeakTS>{TSCermin} && Cer_Quartz_AvgPeakTS<{TSCermax}").Histo1D((
                "hist_FERS_DR_DRSPeak_Cer_Quartz_corr",
                f"Combined energy corrected by DRS Cer Quartz signals;Combined Energy;Counts",
                100, config["xmin_total"]["Mix_sci"], config["xmax_total"]["Mix_sci"]),
                "Energy_DR_corr_Cer_Quartz"
                )

            hist_FERS_DR_filter_DRSPeak_Sci = rdf.Filter(f"Sci_AvgPeakTS>{TSmin} && Sci_AvgPeakTS<{TSmax}").Histo1D((
                "hist_FERS_DR_filter_DRSPeak_Sci",
                f"Combined energy with DRS Sci signals;Combined Energy;Counts",
                100, config["xmin_total"]["Mix_sci"], config["xmax_total"]["Mix_sci"]),
                "Energy_DR"
                )
            hist_FERS_DR_DRSPeak_Sci_corr = rdf.Filter(f"Sci_AvgPeakTS>{TSmin} && Sci_AvgPeakTS<{TSmax}").Histo1D((
                "hist_FERS_DR_DRSPeak_Sci_corr",
                f"Combined energy corrected by DRS Sci signals;Combined Energy;Counts",
                100, config["xmin_total"]["Mix_sci"], config["xmax_total"]["Mix_sci"]),
                "Energy_DR_corr_Sci"
                )
            hists_DRSPeakTS_vs_FERS_EnergySum.append(hist_FERS_DR_filter_DRSPeak_Cer_Plastic)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(hist_FERS_DR_DRSPeak_Cer_Plastic_corr)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(hist_FERS_DR_filter_DRSPeak_Cer_Quartz)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(hist_FERS_DR_DRSPeak_Cer_Quartz_corr)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(hist_FERS_DR_filter_DRSPeak_Sci)
            hists_DRSPeakTS_vs_FERS_EnergySum.append(hist_FERS_DR_DRSPeak_Sci_corr)
            
            
        hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_FERS_Cer_VS_Sci_filter_DRSPeak_Cer_Quartz)
        hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_FERS_Cer_VS_Sci_filter_DRSPeak_Cer_Plastic)
        hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_FERS_Cer_VS_Sci_filter_DRSPeak_Sci)
        hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_FERS_Cer_VS_Sci_DRSPeak_Cer_Quartz_corr)
        hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_FERS_Cer_VS_Sci_DRSPeak_Cer_Plastic_corr)
        hists_DRSPeakTS_vs_FERS_EnergySum.append(h2_FERS_Cer_VS_Sci_DRSPeak_Sci_corr)


    return hists_DRSPeakTS_vs_FERS_EnergySum,proj_DRSPeakTS_vs_FERS_EnergySum



def makeDRSPeakTSVSEnergyCorrPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeakTS_VS_Energy"
    infile_name = f"{rootdir}/drspeak_vs_fers_corr.root"
    infile = ROOT.TFile(infile_name, "READ")
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)
        for cat in ["cer", "sci"]:
            hist_names = [
                f"hist_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}_corr",
                f"hist_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}_corr",
                f"hist_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}_corr",
                ]
            proj_names = [
                f"proj_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}",
                f"proj_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}",
                f"proj_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}",
                ]
            proj_corr_names = [
                f"proj_DRSPeakTS_Cer_Quartz_VS_FERS_{cat}_{gain}_corr",
                f"proj_DRSPeakTS_Cer_Plastic_VS_FERS_{cat}_{gain}_corr",
                f"proj_DRSPeakTS_Sci_VS_FERS_{cat}_{gain}_corr",
                ]
            xtitles = [
                "Cer Quartz peak TS",
                "Cer Plastic peak TS",
                "Sci peak TS"
                ]
            ytitle = f"FERS Energy {cat} {gain}"
            for hist_name,proj_name,proj_corr_name, xtitle in zip(hist_names,proj_names, proj_corr_names, xtitles):
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

                proj = infile.Get(proj_name)
                proj_corr = infile.Get(proj_corr_name)
                extraToDraw = ROOT.TPaveText(0.20, 0.6, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(f"Mean = {proj.GetMean():.3f}, Std = {proj.GetRMS():.3f}")
                extraToDraw.AddText(f"std/mean = {proj.GetRMS()/proj.GetMean():.3f}")
                extraToDraw.AddText(f"Mean corr = {proj_corr.GetMean():.3f}, Std corr = {proj_corr.GetRMS():.3f}")
                extraToDraw.AddText(f"std/mean corr = {proj_corr.GetRMS()/proj_corr.GetMean():.3f}")
                proj.Rebin(4)
                proj_corr.Rebin(4)
                DrawHistos([proj,proj_corr], ["Raw","Corr"], config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"], ytitle, 0, None, "Counts",
                       proj_name,
                       dology=False, drawoptions="HIST", addOverflow=False, addUnderflow=False,
                       outdir=outdir_plots, run_number=run_number, extraToDraw=extraToDraw)
                plots.append(proj_name + ".png")

        hist_names = [
                f"hist_FERS_{gain}_Cer_VS_Sci_filter_DRSPeak_Cer_Quartz",
                f"hist_FERS_{gain}_Cer_VS_Sci_filter_DRSPeak_Cer_Plastic",
                f"hist_FERS_{gain}_Cer_VS_Sci_filter_DRSPeak_Sci",
                f"hist_FERS_{gain}_Cer_VS_Sci_DRSPeak_Cer_Quartz_corr",
                f"hist_FERS_{gain}_Cer_VS_Sci_DRSPeak_Cer_Plastic_corr",
                f"hist_FERS_{gain}_Cer_VS_Sci_DRSPeak_Sci_corr"
                ]
        xtitle = f"FERS Energy Sci {gain}"
        ytitle = f"FERS Energy Cer {gain}"
        for hist_name in hist_names:
            output_name = hist_name.replace("hist_","")
            hist = infile.Get(hist_name)
            extraToDraw = ROOT.TPaveText(0.20, 0.85, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(f"correlation = {hist.GetCorrelationFactor():.3f}")
            DrawHistos([hist], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"], xtitle, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"], ytitle,
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
                       outdir=outdir_plots, addOverflow=False, run_number=run_number,extraToDraw=[extraToDraw])
            plots.append(output_name + ".png")
    

    hist_names = [
        "hist_FERS_DR_filter_DRSPeak_Cer_Quartz",
        "hist_FERS_DR_filter_DRSPeak_Cer_Plastic",
        "hist_FERS_DR_filter_DRSPeak_Sci"
        ]
    hist_names_corr = [
        "hist_FERS_DR_DRSPeak_Cer_Quartz_corr",
        "hist_FERS_DR_DRSPeak_Cer_Plastic_corr",
        "hist_FERS_DR_DRSPeak_Sci_corr"
        ]

    for hist_name, hist_name_corr in zip(hist_names,hist_names_corr):
        hist = infile.Get(hist_name)
        hist_corr = infile.Get(hist_name_corr)
        hist_corr.SetLineColor(ROOT.kRed)
        hist_corr.SetMarkerColor(ROOT.kRed)
        extraToDraw = ROOT.TPaveText(0.20, 0.6, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"Mean = {hist.GetMean():.3f}, Std = {hist.GetRMS():.3f}")
        extraToDraw.AddText(f"std/mean = {hist.GetRMS()/hist.GetMean():.3f}")
        extraToDraw.AddText(f"Mean corr = {hist_corr.GetMean():.3f}, Std corr = {hist_corr.GetRMS():.3f}")
        extraToDraw.AddText(f"std/mean corr = {hist_corr.GetRMS()/hist_corr.GetMean():.3f}")
        hist.Rebin(4)
        hist_corr.Rebin(4)
        xtitle = "FERS Energy Combined"
        output_name = hist_name_corr.replace("hist_","")
        DrawHistos([hist,hist_corr], ["Raw","Corr"], config["xmin_total"]["Mix_sci"], config["xmax_total"]["Mix_sci"], xtitle, 0, None, "Counts",
            output_name,dology=False, drawoptions="HIST", addOverflow=False, addUnderflow=False,
            outdir=outdir_plots, run_number=run_number, extraToDraw=extraToDraw
        )
        plots.append(output_name + ".png")


    output_html = f"{htmldir}/FERSvsDRS//EnergySum_VS_DRSTS_corr.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html



def main():
    makeHists = True
    makePlots = True

    if makeHists:
        global rdf
        rdf = getDRSPeak(rdf, DRSBoards, 450, 550)
        rdf = calibrateDRSPeakTS(rdf, run_number, DRSBoards,
                                 TSminDRS=450, TSmaxDRS=550, threshold=100.0)

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

        hists_DRSTS_Energy,projs_DRSTS_Energy = checkDRSPeakTSVSEnergyCorr(rdf)
        outfile_DRSPeakTS_FERS = ROOT.TFile(
            f"{rootdir}/drspeak_vs_fers_corr.root", "RECREATE")
        for hist in hists_DRSTS_Energy:
            hist.SetDirectory(outfile_DRSPeakTS_FERS)
            hist.Write()
        for proj in projs_DRSTS_Energy:
            proj.SetDirectory(outfile_DRSPeakTS_FERS)
            proj.Write()

    if makePlots:
        output_html_DRSPeakTSVSEnergy = makeDRSPeakTSVSEnergyCorrPlots()
        print(f"DRS Peak TS VS energy plots saved to {output_html_DRSPeakTSVSEnergy}")


if __name__ == "__main__":
    main()

