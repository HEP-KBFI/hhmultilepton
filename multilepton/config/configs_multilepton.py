# coding: utf-8

"""
Configuration of the HH → multi-leptons analysis.
"""

from __future__ import annotations

import os
import re
import itertools
import functools
import yaml
import law
import json
import correctionlib

import order as od

from collections import defaultdict
from scinum import Number

from columnflow.tasks.external import ExternalFile as Ext
from columnflow.util import DotDict, dev_sandbox
from columnflow.columnar_util import ColumnCollection, skip_column
from columnflow.production.cms.top_pt_weight import TopPtWeightConfig
from columnflow.production.cms.dy import DrellYanConfig
from columnflow.production.cms.btag import BTagSFConfig
from columnflow.production.cms.jet import JetIdConfig
from columnflow.production.cms.electron import ElectronSFConfig
from columnflow.production.cms.muon import MuonSFConfig
from columnflow.calibration.cms.tau import TECConfig
from columnflow.calibration.cms.egamma import EGammaCorrectionConfig
from columnflow.config_util import get_root_processes_from_campaign, get_shifts_from_sources
from columnflow.config_util import add_shift_aliases, verify_config_processes 

from multilepton.config.styles import stylize_processes, setup_plot_styles
from multilepton.config.categories import add_categories
from multilepton.config.variables import add_variables
from multilepton.config.met_filters import add_met_filters
from multilepton.config.triggers import add_triggers


logger = law.logger.get_logger(__name__)

# Load analysis configuration
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis.yaml"), "r") as f:
    analysis_data = yaml.load(f, yaml.Loader)


class AnalysisConfig:
    """Helper class to manage analysis configuration from YAML."""
    
    def __init__(self, data):
        self.data = DotDict.wrap(data)
    
    def get_era_key(self, campaign):
        """Get era key for b-tag WPs, luminosity, etc."""
        year = campaign.x.year
        postfix = campaign.x.postfix
        
        if year in [2016, 2017, 2018]: # Run2
            return f"{year}{postfix}" if year == 2016 and postfix == "APV" else str(year)
        else:  # Run 3
            era_map = {
                (2022, ""): "2022",
                (2022, "EE"): "2022EE", 
                (2023, ""): "2023",
                (2023, "BPix"): "2023BPix",
            }
            return era_map.get((year, postfix), str(year))
    
    def get_era(self, campaign):
        year = campaign.x.year
        if year == 2016:
            era= "preVFP" if campaign.has_tag("preVFP") else "postVFP"
        elif year == 2022:
            era = "preEE" if campaign.has_tag("preEE") else "postEE"
        elif year == 2023:
            era = "preBPix" if campaign.has_tag("preBPix") else "postBPix"
        else:
            era = ""
        return f"{year}{era}"

    def get_luminosity(self, campaign):
        """Get luminosity for given campaign."""
        year = campaign.x.year
        lumis = self.data.luminosity.get(str(year)).get("luminosity")
        if isinstance(lumis, dict):
            return lumis.get(get_era(self, campaign))
        else:
            return lumis
   
    def get_process_list(self, process_type="all"):
        """
        Return a flattened list of all 'cmsdb' entries under 'signal' and/or 'background'.
        Args:
            process_type (str): "signal", "background", or "all"
        Returns:
            list[str]: List of process names (from cmsdb entries)
        """
        datasets = self.data.get("datasets", {})
        process_names = []
        categories = ["signal", "background"] if process_type == "all" else [process_type]
    
        for category in categories:
            category_data = datasets.get(category, {})
            if not isinstance(category_data, dict):
                continue
            # Recursively walk through all nested dicts
            def extract_cmsdb_entries(node):
                if isinstance(node, dict):
                    for key, value in node.items():
                        if key == "cmsdb" and isinstance(value, list):
                            process_names.extend(value)
                        else:
                            extract_cmsdb_entries(value)
                elif isinstance(node, list):
                    # In case there are lists of dicts
                    for item in node:
                        extract_cmsdb_entries(item)
            extract_cmsdb_entries(category_data)
        return process_names

    

# Initialize config helper
analysis_cfg = AnalysisConfig(analysis_data)


def pogEraFormat(era):
    """Format era for POG file paths."""
    if any(x in era for x in ['2022', '2023']):
        return era[:4] + '_Summer' + era.replace('20', '')
    else:
        return era.replace("UL", "") + "_UL"


def localizePOGSF(era, POG, fileName):
    """Localize POG scale factor files."""
    subdir = pogEraFormat(str(era))
    return os.path.join("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration", "POG", POG, subdir, fileName)


def nested_dict():
    """Recursively create nested defaultdicts."""
    return defaultdict(nested_dict)


#https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
def bTagWorkingPoints(year, run, campaign):
    fileName = localizePOGSF(year, "BTV", "btagging.json.gz")
    logger.info(f'... getting working points and discr cuts from : {fileName}')
    ceval = correctionlib.CorrectionSet.from_file(fileName)
    btagging = nested_dict()
    if run == 2:
        taggers = ["deepJet", "deepcsv", "particleNetMD"]
        eras = ["2016APV", "2016", "2017", "2018"]
    elif run == 3:
        taggers = ["deepJet", "particleNet", "robustParticleTransformer", "particleNetMD"]
        eras = ["2022", "2022EE", "2023", "2023BPix"]
    else:
        raise ValueError(f"Unsupported run: {run}")

    era = f"{year}{campaign.x.postfix}"
    for tagger in taggers:
        for wp in ['L', 'M', 'T', 'XT', 'XXT']:
            try:
                btagging[tagger][wp][era] = ceval[f"{tagger.replace('MD', '')}_wp_values"].evaluate(wp)
            except Exception as e:
                logger.warning(f"Failed to evaluate {tagger} {wp} for {era}: {e}")
    # Optionally convert defaultdicts to normal dicts for output
    return json.loads(json.dumps(btagging))


def build_leaf_processes(procs, base_name, jet_configs):
    """Build leaf processes for a given base process name."""
    leaf_processes = []
    for n_jets, pt_bins in jet_configs:
        if n_jets == 0:
            leaf_processes.append(procs.get(f"{base_name}_0j"))
        elif pt_bins:
            leaf_processes.extend(
                procs.get(f"{base_name}_{n_jets}j_pt{pt}")
                for pt in pt_bins
            )
    # Add >=3 jets process
    leaf_processes.append(procs.get(f"{base_name}_ge3j"))
    return leaf_processes


def convert_dataset_to_process(dataset, campaign, all_processes_from_campaign):
    process = dataset
    for production in ['_powheg', '_amcatnlo', '_pythia', '_madgraph']:
        if production in dataset:
            process = dataset.replace(production, '')
    # Find matching process and return its id
    id = None
    for proc in all_processes_from_campaign:
        if process == proc.name:
            id = proc.id
            break  # <-- exit the loop immediately when found
    if id is None:
        #raise ValueError(f"Will skip ... No matching process '{process}' found for dataset in campaign '{campaign.name}'")
        logger.warning(f"Will skip ... No matching process '{process}' found for dataset in campaign '{campaign.name}'")
    return process, id


def add_config(
    analysis: od.Analysis,
    campaign: od.Campaign,
    config_name: str | None = None,
    config_id: int | None = None,
    limit_dataset_files: int | None = None,
) -> od.Config:
    
    # gather campaign data
    run = campaign.x.run
    year = campaign.x.year
 
    # --- basic configuration validations ---
    if run not in {2, 3}:
        raise ValueError(f"Invalid run: {run}. Expected 2 or 3.")
    
    valid_years = {2016, 2017, 2018, 2022, 2023, 2024, 2025}
    if year not in valid_years:
        raise ValueError(f"Invalid year: {year}. Must be one of {sorted(valid_years)}.")
    
    # get all root processes
    all_processes_from_campaign = get_root_processes_from_campaign(campaign)
    for proc in all_processes_from_campaign:
        print( proc)#, proc.name, proc.id) 
    
    # create a config by passing the campaign
    cfg = od.Config(
        name=config_name,
        id=config_id,
        campaign=campaign,
    )
    
    #=============================================
    # helpers
    #=============================================
    def set_luminosity(campaign, year, analysis_data):
        year_data = analysis_data["years"].get(year)
        if not year_data:
            raise ValueError(f"Year {year} not found in analysis.yaml")
        # detect campaign tag (e.g. preVFP, postEE, etc.)
        tag = next((t for t in ["preVFP", "postVFP", "preEE", "postEE", "preBPix", "postBPix"]
                    if campaign.has_tag(t)), None)
        lumi_info = year_data["luminosity"]
        if isinstance(lumi_info, list):
            lumi_map = {list(d.keys())[0]: list(d.values())[0] for d in lumi_info}
            key = f"{year}{tag}" if tag else list(lumi_map.keys())[0]
            lumi_value = lumi_map.get(key)
        else:
            lumi_value = lumi_info
        lumi_unc = year_data.get("luminosity-uncertainties", {})
        return lumi_value, lumi_unc 
   
    def _names_from_tag(tag):
        return [s.name for s in cfg.shifts if s.has_tag(tag)]

    def ConfigureMuons(cfg, run, year, campaign):
        if run == 2:
            cfg.x.muon_sf_names = MuonSFConfig(
                correction="NUM_TightRelIso_DEN_TightIDandIPCut"
            )
        elif run == 3:
            cfg.x.muon_sf_names = MuonSFConfig(
                correction="NUM_TightPFIso_DEN_TightID"
            )
            cfg.x.muon_trigger_sf_names = MuonSFConfig("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight")
            cfg.x.single_trigger_muon_data_effs_cfg = MuonSFConfig("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_DATAeff")
            cfg.x.single_trigger_muon_mc_effs_cfg = MuonSFConfig("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_MCeff")
            cfg.x.cross_trigger_muon_data_effs_cfg = MuonSFConfig("NUM_IsoMu20_DEN_CutBasedIdTight_and_PFIsoTight_DATAeff")
            cfg.x.cross_trigger_muon_mc_effs_cfg = MuonSFConfig("NUM_IsoMu20_DEN_CutBasedIdTight_and_PFIsoTight_MCeff")
        return cfg 

    def ConfigureElectrons(cfg, run, year, campaign):
        """ Run 2: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaULTagAndProbe
            Run 3: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRun3Recommendations
        """
        EGMcorrection = {
                "2016APV": "preVFP",
                "2016": "postVFP",
                "2017": "",
                "2018": "",
                "2022EE": "Re-recoE+PromptFG", 
                "2022": "Re-recoBCD",
                "2023BPix": "PromptD",
                "2023": "PromptC"
                }
        e_postfix = EGMcorrection.get(f"{year}{campaign.x.postfix}")
        e_prefix = 'UL-' if run == 2 else ''
        cfg.x.electron_sf_names = ElectronSFConfig(
            correction="{e_prefix}Electron-ID-SF",
            campaign=f"{year}{e_postfix}",
            working_point="wp80iso",
        )
        cfg.x.eec = EGammaCorrectionConfig(
            correction_set="Scale",
            value_type="total_correction",
            uncertainty_type="total_uncertainty",
        )
        cfg.x.eer = EGammaCorrectionConfig(
            correction_set="Smearing",
            compound=False,
            value_type="rho",
            uncertainty_type="err_rho",
        )
        # Define HLT paths for easier maintenance
        hlt_single, hlt_cross = "HLT_SF_Ele30_TightID", "HLT_SF_Ele24_TightID"

        # Common helper for trigger configs
        def make_el_trigger_cfg(corr, path, campaign=f"{year}{e_postfix}", suffix=""):
            return ElectronSFConfig(correction=f"Electron-HLT-{corr}{suffix}",
                                    campaign=campaign, hlt_path=path)

        cfg.x.electron_trigger_sf_names = make_el_trigger_cfg("SF", hlt_single)
        cfg.x.single_trigger_electron_data_effs_cfg = make_el_trigger_cfg("DataEff", hlt_single)
        cfg.x.single_trigger_electron_mc_effs_cfg = make_el_trigger_cfg("McEff", hlt_single)
        cfg.x.cross_trigger_electron_data_effs_cfg = make_el_trigger_cfg("DataEff", hlt_cross)
        cfg.x.cross_trigger_electron_mc_effs_cfg = make_el_trigger_cfg("McEff", hlt_cross)

        # --- Electron Energy Corrections (EEC/EER) ----------------------------------------------
        if year == 2022:
            e_tag = {"": "preEE", "EE": "postEE"}[campaign.x.postfix]
        elif year == 2023:
            e_tag = {"": "preBPIX", "BPix": "postBPIX"}[campaign.x.postfix]
        
        cfg.x.eec = EGammaCorrectionConfig(
            correction_set=f"EGMScale_Compound_Ele_{year}{e_tag}",
            value_type="scale",
            uncertainty_type="escale",
            compound=True,
        )
        cfg.x.eer = EGammaCorrectionConfig(
            correction_set=f"EGMSmearAndSyst_ElePTsplit_{year}{e_tag}",
            value_type="smear",
            uncertainty_type="esmear",
        )
        return cfg

    def ConfigureTaus(cfg, run, campaign):
        """
        Configure tau ID, TEC (Tau Energy Calibration), and trigger settings.
        """
        tau_taggers = {
            2: "DeepTau2017v2p1",
            3: "DeepTau2018v2p5",
        }
    
        cfg.x.tau_tagger = tau_taggers.get(run)
        corrector_kwargs = {"wp": "Medium", "wp_VSe": "VVLoose"} if run == 3 else {}
        cfg.x.tec = TECConfig(tagger=cfg.x.tau_tagger, corrector_kwargs=corrector_kwargs)
    
        # --- Tau ID working points
        # Legacy (campaign.x.version < 10) vs New format (>=10)
        if campaign.x.version < 10:
            wp_values_mu = {"vloose": 1, "loose": 2, "medium": 4, "tight": 8}
            wp_values_jet_or_e = {
                "vvvloose": 1, "vvloose": 2, "vloose": 4,
                "loose": 8, "medium": 16, "tight": 32,
                "vtight": 64, "vvtight": 128,
            }
        else:
            wp_values_mu = {"vloose": 1, "loose": 2, "medium": 3, "tight": 4}
            wp_values_jet_or_e = {
                "vvvloose": 1, "vvloose": 2, "vloose": 3,
                "loose": 4, "medium": 5, "tight": 6,
                "vtight": 7, "vvtight": 8,
            }
    
        cfg.x.tau_id_working_points = DotDict.wrap({
            "tau_vs_e": wp_values_jet_or_e,
            "tau_vs_jet": wp_values_jet_or_e,
            "tau_vs_mu": wp_values_mu,
        })
    
        # --- Tau trigger working points
        cfg.x.tau_trigger_working_points = DotDict.wrap({
            "id_vs_jet_v0": "VVLoose",
            "id_vs_jet_gv0": ("Loose", "VVLoose"),
            "id_vs_mu_single": "Tight",
            "id_vs_mu_cross": "VLoose",
            "id_vs_e_single": "VVLoose",
            "id_vs_e_cross": "VVLoose",
            "trigger_corr": "VVLoose",
        })
        # --- Tau trigger correctors
        cfg.x.tau_trigger_corrector = "tau_trigger"
        cfg.x.tau_trigger_corrector_cclub = "tauTriggerSF"
        return cfg

    def ConfigureJets(cfg, year, run, campaign):
        """
        Configure JEC and JER campaigns and versions for Run 2 and Run 3.
        References:
            - Run 2: https://cms-jerc.web.cern.ch/Recommendations/#run-2
                      https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=204
                      https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=109
            - Run 3: https://cms-jerc.web.cern.ch/Recommendations/#2022
        """
        newyear = year%100
        if run == 2:
            jec_version_map = {2016: "V7", 2017: "V5", 2018: "V5"}
            jer_version_map = {2016: "V3", 2017: "V2", 2018: "V2"}
            jec_campaign = f"Summer19UL{newyear}{campaign.x.postfix}"
            jer_campaign = f"Summer{'20' if year == 2016 else '19'}UL{newyear}{campaign.x.postfix}"
            return {
                "jec_campaign": jec_campaign,
                "jec_version": jec_version_map[year],
                "jer_campaign": jer_campaign,
                "jer_version": "JR" + jer_version_map[year],
                "jet_type": "AK4PFchs",
                "data_per_era": False,
            }
        elif run == 3:
            jerc_postfix = {
                (2022, ""): "_22Sep2023",
                (2022, "EE"): "_22Sep2023",
                (2023, ""): "Prompt23",
                (2023, "BPix"): "Prompt23",
            }.get((year, campaign.x.postfix))
            jec_version_map = {
                (2022, ""): "V2",
                (2022, "EE"): "V2",
                (2023, ""): "V2",
                (2023, "BPix"): "V3",
            }
            if not jerc_postfix:
                raise ValueError(f"Unsupported JERC configuration for Run 3: year={year}, postfix={campaign.x.postfix}")
            jec_campaign = f"Summer{newyear}{campaign.x.postfix}{jerc_postfix}"
            jer_campaign = f"Summer{newyear}{campaign.x.postfix}{jerc_postfix}"
            # Add special Run fragment for 2023
            if year == 2023:
                jer_campaign += f"_Run{'Cv1234' if campaign.has_tag('preBPix') else 'D'}"
            return {
                "jec_campaign": jec_campaign,
                "jec_version": jec_version_map[(year, campaign.x.postfix)],
                "jer_campaign": jer_campaign,
                "jer_version": "JR" + {2022: "V1", 2023: "V1"}[year],
                "jet_type": "AK4PFPuppi",
                "data_per_era": year == 2022,  # 2022 JEC depends on era
            }

    def ConfigureLFNS(cfg, limit_dataset_files=None):
        """
        Configure custom methods for retrieving dataset LFNs depending on campaign settings.
        """
        cfg.x.get_dataset_lfns = None
        cfg.x.get_dataset_lfns_sandbox = None
        
        # Handle special campaign type: "custom" with "creator" == "uhh"
        campaign_custom = cfg.campaign.x("custom", {})
        if campaign_custom.get("creator") != "uhh":
            return cfg  # No custom configuration needed
    
        def get_multileptons_dataset_lfns(dataset_inst: od.Dataset, shift_inst: od.Shift, dataset_key: str) -> list[str]:
            """
            Retrieve LFNs for a given dataset under the UHH custom campaign convention.
            """
            try:
                _, dataset_id, full_campaign, tier = dataset_key.split("/")
                main_campaign, sub_campaign = full_campaign.split("-", 1)
            except ValueError:
                raise ValueError(f"Invalid dataset key format: {dataset_key}")
    
            path = f"store/{dataset_inst.data_source}/{main_campaign}/{dataset_id}/{tier}/{sub_campaign}/0"
            # Determine filesystem and directory target class
            custom_name = campaign_custom.get("name")
            remote_fs = f"wlcg_fs_{custom_name}"
            local_fs = f"local_fs_{custom_name}"
            dir_cls = law.wlcg.WLCGDirectoryTarget
            fs_to_use = remote_fs
    
            if law.config.has_section(local_fs):
                base = law.target.file.remove_scheme(law.config.get_expanded(local_fs, "base"))
                if os.path.exists(base):
                    dir_cls = law.LocalDirectoryTarget
                    fs_to_use = local_fs
    
            lfn_base = dir_cls(path, fs=fs_to_use)
            # Retrieve all ROOT files and convert them to LFNs
            lfns = [
                "/" + lfn_base.child(fname, type="f").path.lstrip("/")
                for fname in lfn_base.listdir(pattern="*.root")
            ]
            return sorted(lfns)
        print( get_multileptons_dataset_lfns )
        # Attach the retrieval method and related configuration
        cfg.x.get_dataset_lfns = get_multileptons_dataset_lfns
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
        cfg.x.get_dataset_lfns_remote_fs = lambda dataset_inst: [
            f"local_fs_{campaign_custom['name']}",
            f"wlcg_fs_{campaign_custom['name']}",
        ]
        return cfg

    def add_external(name, value):
        if isinstance(value, dict):
            value = DotDict.wrap(value)
        cfg.x.external_files[name] = value

    def register_shift_pair(cfg, base_name, base_id, aliases=None, tags=None, aux=None, step=1):
        """Register up/down shifts with optional aliases, tags, and aux data."""
        cfg.add_shift(name=f"{base_name}_up", id=base_id, type="shape", tags=tags or set(), aux=aux)
        cfg.add_shift(name=f"{base_name}_down", id=base_id + step, type="shape", tags=tags or set(), aux=aux)
        if aliases:
            add_shift_aliases(cfg, base_name, aliases)
    
    def find_match_era(**kwargs):
        """Helper to enable processes/datasets only for specific era."""
        return (
            (kwargs.get('run') is None or campaign.x.run in law.util.make_set(kwargs.get('run'))) and
            (kwargs.get('year') is None or campaign.x.year in law.util.make_set(kwargs.get('year'))) and
            (kwargs.get('postfix') is None or campaign.x.postfix in law.util.make_set(kwargs.get('postfix'))) and
            (kwargs.get('tag') is None or campaign.has_tag(kwargs.get('tag'), mode=any)) and
            (kwargs.get('nano') is None or campaign.x.version in law.util.make_set(kwargs.get('nano')))
        )

    def if_era(values=None, **kwargs):
        return list(filter(bool, values or [])) if find_match_era(**kwargs) else []

    def if_not_era(values=None, **kwargs):
        return list(filter(bool, values or [])) if not find_match_era(**kwargs) else []

    def if_in_config_id(ids, values):
        return list(filter(bool, values or [])) if cfg.id in ids else []

    def if_not_in_config_id(ids, values):
        return list(filter(bool, values or [])) if cfg.id not in ids else [] 
    
    #=============================================
    # configure some default objects
    #=============================================
    cfg.x.default_calibrator = "default"
    cfg.x.default_selector = "default"
    cfg.x.default_reducer = "default"
    cfg.x.default_producer = "default"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "default_no_shifts"
    cfg.x.default_categories = ("all",)
    cfg.x.default_variables = ("njet", "nbtag", "res_pdnn_hh", "res_dnn_hh")
    cfg.x.default_hist_producer = "default"
    cfg.x.external_files = DotDict() 
    
    # get btag working points 
    cfg.x.btag_working_points = bTagWorkingPoints(year, run, campaign)
    
    btagJECsources = analysis_data.get("btag_sf_jec_sources", [])
    btagJECsources += [f"Absolute_{year}", f"BBEC1_{year}", f"EC2_{year}", f"HF_{year}", f"RelativeSample_{year}", ""] 
    cfg.x.btag_sf_jec_sources = btagJECsources
    
    lumi_value, lumi_unc = set_luminosity(campaign, year, analysis_data)
    #cfg.x.luminosity = Number(lumi_value, lumi_unc)
   
    ConfigureJets(cfg, year, run, campaign)
    ConfigureTaus(cfg, run, campaign)
    ConfigureElectrons(cfg, run, year, campaign)
    ConfigureMuons(cfg, run, year, campaign)
    ConfigureLFNS(cfg, limit_dataset_files)
    #=============================================
    # LFNS, processes and datasets - using YAML configuration
    #=============================================
    dataset_names = []
    datasets_config = analysis_data.get("datasets", {})
    process_names = analysis_cfg.get_process_list("all")
   
    # Loop over signal and background
    for dtype in ['signal', 'background']:
        for group, dCfg in datasets_config.get(dtype, {}).items():
            # Handle special case (QCD)
            if group == "qcd":
                from cmsdb.processes.qcd import qcd
                cfg.add_process(qcd)
                dataset_names.extend(qcd.names if hasattr(qcd, "names") else [])
                continue
            # Recursively find cmsdb lists
            def extract_cmsdb_entries(node):
                if isinstance(node, dict):
                    if "cmsdb" in node:
                        return node["cmsdb"]
                    result = []
                    for value in node.values():
                        result.extend(extract_cmsdb_entries(value))
                    return result
                return []
            datasets = extract_cmsdb_entries(dCfg)
            tags = dCfg.get("tag", [])
             
            for dataset_name in datasets:
                dataset_names.append(dataset_name)
                proc, id = convert_dataset_to_process(dataset_name, campaign, all_processes_from_campaign)
                if id is None or not campaign.has_dataset(dataset_name):
                    continue
                cfg.add_process(proc, id)
                dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
                # Add tags to the process, not the string
                for tag in (t for t in law.util.make_set(tags) if t is not None):
                    dataset.add_tag(tag)
                if limit_dataset_files:
                    for info in dataset.info.values():
                        info.n_files = min(info.n_files, limit_dataset_files)
    
    # Add data
    streams = datasets_config["data"]["streams"]
    for year, year_cfg in datasets_config["data"].items():
        if year == "streams":
            continue
    
        # year_cfg is a list of tag blocks
        for tag_block in year_cfg:
            # tag_block is a dict like {"preEE": {...}}
            for tag, tag_cfg in tag_block.items():
                periods = tag_cfg["periods"]
    
                # Example: add dataset names
                dataset_names += [
                    *if_not_in_config_id(
                        ids=[5012, 6012, 7012, 8012],
                        values=[
                            *if_era(
                                year=year,
                                **({"tag": tag} if tag else {}),  # only add tag if it exists
                                values=[f"data_{stream}_{period}" for stream in streams for period in periods])]
                            )]
                # Optional: special tag for broken MET filter in 2022
                # bad ecalBadCalibFilter MET filter in 2022 data
                # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=172#ECal_BadCalibration_Filter_Flag
                # https://cms-talk.web.cern.ch/t/noise-met-filters-in-run-3/63346/5
                if year == 2022 and dataset.is_data and dataset.x.era in "FG":
                    dataset.add_tag("broken_ecalBadCalibFilter")
                if limit_dataset_files:
                    for info in dataset.info.values():
                        info.n_files = min(info.n_files, limit_dataset_files)
    
    # Configure colors, labels, etc
    stylize_processes(cfg)
   
    #=============================================
    # Jet Energy Corrections (JEC) and Jet Energy Resolution (JER)
    #=============================================
    jecjerdb = ConfigureJets(cfg, year, run, campaign)
    cfg.x.jec = DotDict.wrap({
        "Jet": {
            "campaign": jecjerdb["jec_campaign"],
            "version": jecjerdb["jec_version"],
            "data_per_era": jecjerdb["data_per_era"],
            "jet_type": jecjerdb["jet_type"],
            "levels": ["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],
            "levels_for_type1_met": ["L1FastJet"],
            "uncertainty_sources": analysis_data['jec_sources'],
        },
    })
    
    cfg.x.jer = DotDict.wrap({
        "Jet": {
            "campaign": jecjerdb["jer_campaign"],
            "version": jecjerdb["jer_version"],
            "jet_type": jecjerdb["jet_type"],
        },
    })
    
    cfg.x.jet_id = JetIdConfig( 
        corrections={
            "AK4PUPPI_Tight": 2,
            "AK4PUPPI_TightLeptonVeto": 3,
        })

    cfg.x.fatjet_id = JetIdConfig(
        corrections={
            "AK8PUPPI_Tight": 2,
            "AK8PUPPI_TightLeptonVeto": 3,
        })
    
    cfg.x.jet_trigger_corrector = "jetlegSFs" 
    
    #=============================================
    # met settings
    #=============================================
    if run == 2:
        cfg.x.met_name = "MET"
        cfg.x.raw_met_name = "RawMET"
        cfg.x.met_phi_correction_set = r"{variable}_metphicorr_pfmet_{data_source}"
    elif run == 3:
        cfg.x.met_name = "PuppiMET"
        cfg.x.raw_met_name = "RawPuppiMET"

    #=============================================
    # b-tag working points 
    #=============================================
    cfg.x.btag_sf_deepjet = BTagSFConfig(
        correction_set="deepJet_shape",
        jec_sources=cfg.x.btag_sf_jec_sources,
        discriminator="btagDeepFlavB",
    )
    if run == 3:
        cfg.x.btag_sf_pnet = BTagSFConfig(
            correction_set="particleNet_shape",
            jec_sources=cfg.x.btag_sf_jec_sources,
            discriminator="btagPNetB",
        )

    #=============================================
    # top pt reweighting
    # https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting?rev=31
    #=============================================
    cfg.x.top_pt_weight = TopPtWeightConfig(
        params={
            "a": 0.0615,
            "a_up": 0.0615 * 1.5,
            "a_down": 0.0615 * 0.5,
            "b": -0.0005,
            "b_up": -0.0005 * 1.5,
            "b_down": -0.0005 * 0.5,
        },
        pt_max=500.0,
    )

    #=============================================
    # dy reweighting and recoil 
    #=============================================
    if run == 3:
        era = analysis_cfg.get_era(campaign)
        
        # dy reweighting
        # https://cms-higgs-leprare.docs.cern.ch/htt-common/DY_reweight
        cfg.x.dy_weight_config = DrellYanConfig(
            era=era,
            order="NLO",
            correction="DY_pTll_reweighting",
            unc_correction="DY_pTll_reweighting_N_uncertainty",
        )

        # dy boson recoil correction
        # https://cms-higgs-leprare.docs.cern.ch/htt-common/V_recoil
        cfg.x.dy_recoil_config = DrellYanConfig(
            era=era,
            order="NLO",
            correction="Recoil_correction_Rescaling",
            unc_correction="Recoil_correction_Uncertainty",
        ) 

    #=============================================
    # shifts
    #=============================================
    # load JEC sources
    all_jec_sources = analysis_data.get("jec_sources", [])
    
    # nominal + simple shifts
    simple_shifts = [
        # (name, base_id, aliases, tags, aux)
        #("nominal", 0, None, None),
        ("tune", 1, None, {"disjoint_from_nominal"}, None),
        ("hdamp", 3, None, {"disjoint_from_nominal"}, None),
        ("mtop", 5, None, {"disjoint_from_nominal"}, None),
        ("minbias_xs", 7, None, None, {
            "pu_weight": "pu_weight_{name}",
            "normalized_pu_weight": "normalized_pu_weight_{name}",
        }),
        ("top_pt", 9, None, None, {"top_pt_weight": "top_pt_weight_{direction}"}),
    ]
    
    for name, base_id, aliases, tags, aux in simple_shifts:
        register_shift_pair(cfg, name, base_id, aliases, tags, aux)
    
    # JEC sources
    for jec_source in cfg.x.jec.Jet.uncertainty_sources:
        idx = all_jec_sources.index(jec_source)
        jec_id = 5000 + 2 * idx
        jec_aliases = {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
            f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
        }
        register_shift_pair(cfg, f"jec_{jec_source}", jec_id, jec_aliases, {"jec"}, {"jec_source": jec_source})
    
        # link btag-related JEC sources
        if ("" if jec_source == "Total" else jec_source) in cfg.x.btag_sf_jec_sources:
            btag_aliases = {
                "normalized_btag_deepjet_weight": "normalized_btag_deepjet_weight_{name}",
                "normalized_njet_btag_deepjet_deepjet_weight": "normalized_njet_btag_deepjet_weight_{name}",
                "normalized_btag_pnet_weight": "normalized_btag_pnet_weight_{name}",
                "normalized_njet_btag_pnet_weight": "normalized_njet_btag_pnet_weight_{name}",
            }
            add_shift_aliases(cfg, f"jec_{jec_source}", btag_aliases)
    
    # JER
    jer_aliases = {
        "Jet.pt": "Jet.pt_{name}",
        "Jet.mass": "Jet.mass_{name}",
        f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
        f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
    }
    register_shift_pair(cfg, "jer", 6000, jer_aliases, {"jer"})
    
    # TEC shifts 
    for i, (match, dm) in enumerate(itertools.product(["jet", "e"], [0, 1, 10, 11])):
        tec_aliases = {
            "Tau.pt": "Tau.pt_{name}",
            "Tau.mass": "Tau.mass_{name}",
            f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
            f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
        }
        register_shift_pair(cfg, f"tec_{match}_dm{dm}", 20 + 2 * i, tec_aliases, {"tec"})
    
    # TAU uncertainties 
    cfg.x.tau_unc_names = [
        "jet_dm0", "jet_dm1", "jet_dm10", "jet_dm11",
        "e_barrel", "e_endcap",
        "mu_0p0To0p4", "mu_0p4To0p8", "mu_0p8To1p2", "mu_1p2To1p7", "mu_1p7To2p3",
    ]
    for i, unc in enumerate(cfg.x.tau_unc_names):
        register_shift_pair(cfg, f"tau_{unc}", 50 + 2 * i, {"tau_weight": f"tau_weight_{unc}_{{direction}}"})
    
    # Electron, muon, and energy corrections
    register_shift_pair(cfg, "e", 90, {"electron_weight": "electron_weight_{direction}"})
    if run == 3 and year == 2022:
        logger.debug("adding ees and eer shifts")
        register_shift_pair(cfg, "ees", 92, {"Electron.pt": "Electron.pt_scale_{direction}"}, {"eec"})
        register_shift_pair(cfg, "eer", 94, {"Electron.pt": "Electron.pt_res_{direction}"}, {"eer"})
    register_shift_pair(cfg, "mu", 100, {"muon_weight": "muon_weight_{direction}"})
    
    # b-tag uncertainties
    cfg.x.btag_unc_names = [
        "hf", "lf",
        f"hfstats1_{year}", f"hfstats2_{year}",
        f"lfstats1_{year}", f"lfstats2_{year}",
        "cferr1", "cferr2",
    ]
    for i, unc in enumerate(cfg.x.btag_unc_names):
        btag_aliases = {
            "normalized_btag_deepjet_weight": f"normalized_btag_deepjet_weight_{unc}_{{direction}}",
            "normalized_njet_btag_deepjet_weight": f"normalized_njet_btag_deepjet_weight_{unc}_{{direction}}",
        }
        register_shift_pair(cfg, f"btag_{unc}", 110 + 2 * i, btag_aliases)
    
    # LHE variations
    lhe_shifts = {
        "pdf": 130,
        "murmuf": 140,
        "isr": 150,
        "fsr": 155,
    }
    for name, base_id in lhe_shifts.items():
        aliases = {
            f"{name}_weight": f"{name}_weight_{{direction}}",
            f"normalized_{name}_weight": f"normalized_{name}_weight_{{direction}}",
        }
        register_shift_pair(cfg, name, base_id, aliases, {"lhe_weight"} if name in ["pdf", "murmuf"] else None)
    
    # trigger scale factors
    trigger_legs = ["e", "mu", "tau_dm0", "tau_dm1", "tau_dm10", "tau_dm11", "jet"]
    for i, leg in enumerate(trigger_legs):
        register_shift_pair(cfg, f"trigger_{leg}", 180 + 2 * i)
   
    #=============================================
    # Add scale-factors from correction lib 
    #=============================================
    if run == 2:
        tauPOGJsonFile = "tau.json.gz"
        metPOGJsonFile = "met.json.gz"
    elif run == 3: # nasty names, workaround 
        if year == 2022:
            met_pog_suffix = f"{year}_{year}{'' if campaign.has_tag('preEE') else 'EE'}"
            tau_pog_suffix = f"{'pre' if campaign.has_tag('preEE') else 'post'}EE"
        elif year == 2023:
            met_pog_suffix = f"_{year}_{year}{'' if campaign.has_tag('preBPix') else 'BPix'}"
            tau_pog_suffix = f"{'pre' if campaign.has_tag('preBPix') else 'post'}BPix"
        tauPOGJsonFile = f"tau_DeepTau2018v2p5_{year}_{tau_pog_suffix}.json.gz"
        metPOGJsonFile = f"met_xyCorrections_{met_pog_suffix}.json.gz"
        triggerFile = analysis_data["external_files"]["trigger"]
        campaign_tag = ""
        for tag in ("preEE", "postEE", "preBPix", "postBPix"):
            if campaign.has_tag(tag, mode=any):
                if campaign_tag:
                    raise ValueError(f"Multiple campaign tags found: {cfg.x.campaign_tag} and {tag}")
                campaign_tag = tag

    # common files
    # (versions in the end are for hashing in cases where file contents changed but paths did not)
    goldenFile = analysis_data['years'][year]["certified_lumi_file"]
    normtagFile = analysis_data['years'][year]["normtag"]

    add_external("lumi", {"golden":(goldenFile, "v1"), "normtag": (normtagFile, "v1")})
    add_external("pu_sf", (localizePOGSF(year, "LUM", "puWeights.json.gz"), "v1"))
    add_external("jet_jerc", (localizePOGSF(year, "JME", "jet_jerc.json.gz"), "v1"))
    add_external("jet_veto_map", (localizePOGSF(year, "JME", "jetvetomaps.json.gz"), "v1"))
    add_external("met_phi_corr", (localizePOGSF(year, "JME", f"{metPOGJsonFile}"), "v1"))
    add_external("btag_sf_corr", (localizePOGSF(year, "BTV", "btagging.json.gz"), "v1"))
    add_external("muon_sf", (localizePOGSF(year, "MUO", "muon_Z.json.gz"), "v1"))
    add_external("electron_sf", (localizePOGSF(year, "EGM", "electron.json.gz"), "v1"))
    add_external("tau_sf", (localizePOGSF(year, "TAU", f"{tauPOGJsonFile}"), "v1"))

    # run specific files
    if run == 2:
        add_external("tau_trigger_sf", (localizePOGSF(year, "TAU", "tau.json.gz"), "v1"))
        # hh-btag repository with TF saved model directories trained on Run2 UL samples
        add_external("hh_btag_repo", Ext(
            "/afs/cern.ch/work/m/mrieger/public/hbt/external_files/hh-btag-master-d7a71eb3.tar.gz",
            subpaths=DotDict(even="hh-btag-master/models/HHbtag_v2_par_0", odd="hh-btag-master/models/HHbtag_v2_par_1"),  # noqa: E501
            version="v2",
        ))
    elif run == 3:
        # electron energy correction and smearing
        add_external("electron_ss", (localizePOGSF(year, "EGM", "electronSS_EtDependent.json.gz"), "v1"))
        # updated jet id
        add_external("jet_id", (localizePOGSF(year, "JME", "jetid.json.gz"), "v1"))
        # hh-btag repository with TF saved model directories trained on 22+23 samples using pnet
        add_external("hh_btag_repo", Ext(
            "/afs/cern.ch/work/m/mrieger/public/hbt/external_files/hh-btag-master-d7a71eb3.tar.gz",
            subpaths=DotDict(even="hh-btag-master/models/HHbtag_v3_par_0", odd="hh-btag-master/models/HHbtag_v3_par_1"),  # noqa: E501
            version="v3",
        ))

        # dy weight and recoil corrections
        add_external("dy_weight_sf", ("/afs/cern.ch/work/m/mrieger/public/mirrors/external_files/DY_pTll_weights_v2.json.gz", "v1"))  # noqa: E501
        add_external("dy_recoil_sf", ("/afs/cern.ch/work/m/mrieger/public/mirrors/external_files/Recoil_corrections_v2.json.gz", "v1"))  # noqa: E501

        # trigger scale factors
        trigger_sf_internal_subpath = "AnalysisCore-59ae66c4a39d3e54afad5733895c33b1fb511c47/data/TriggerScaleFactors"
        add_external("trigger_sf", Ext(
            triggerFile,
            subpaths=DotDict(
                muon=f"{trigger_sf_internal_subpath}/{year}{campaign_tag}/temporary_MuHlt_abseta_pt.json",
                cross_muon=f"{trigger_sf_internal_subpath}/{year}{campaign_tag}/CrossMuTauHlt.json",
                electron=f"{trigger_sf_internal_subpath}/{year}{campaign_tag}/electronHlt.json",
                cross_electron=f"{trigger_sf_internal_subpath}/{year}{campaign_tag}/CrossEleTauHlt.json",
                tau=f"{trigger_sf_internal_subpath}/{year}{campaign_tag}/tau_trigger_DeepTau2018v2p5_{year}{tau_pog_suffix}.json",
                jet=f"{trigger_sf_internal_subpath}/{year}{campaign_tag}/ditaujet_jetleg_SFs_{campaign_tag}.json",
            ),
            version="v1",
        ))

    #=============================================
    # reductions
    #=============================================
    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        # !! note that this set is used by the cf_default reducer
        "cf.ReduceEvents": {
            # mandatory
            ColumnCollection.MANDATORY_COFFEA,
            # event info
            "deterministic_seed",
            # object info
            "Jet.{pt,eta,phi,mass,hadronFlavour,puId,hhbtag,btag*,nConstituents,deterministic_seed}",
            "NonCleanedJet.{pt,eta,phi,mass,hadronFlavour,puId,hhbtag,btag*,nConstituents,deterministic_seed}",
            "VBFJet.{pt,eta,phi,mass,hadronFlavour,puId,hhbtag,btag*,nConstituents,deterministic_seed}",
            "FatJet.*",
            "SubJet{1,2}.*",
            "Electron.*",
            "ElectronLoose.*",
            "ElectronTight.*",
            "Muon.*",
            "MuonLoose.*",
            "MuonTight.*",
            "Tau.*",
            "TauNoID.*",
            "TauIso.*",
            "GenPart*",
            f"{cfg.x.met_name}.{{pt,phi,significance,covXX,covXY,covYY}}",
            "PV.npvs",
            "HLT.*",
            # keep all columns added during selection and reduction, but skip cutflow features
            ColumnCollection.ALL_FROM_SELECTOR,
            skip_column("cutflow.*"),
        },
        "cf.MergeSelectionMasks": {
            "cutflow.*",
        },
        "cf.UniteColumns": {
            "*", *skip_column("*_{up,down}"),
        },
    })
    
    #=============================================
    # Add event weights configuration
    # Each key corresponds to a possible event weight column.
    # Each value is a list of shift dependencies (from systematics sources).
    # This mapping is used by weight producers later in the workflow.
    #=============================================
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    
    # --- Global event weight configuration ---
    base_event_weights = {
        "pdf_weight":                 get_shifts("pdf"),
        "murmuf_weight":              get_shifts("murmuf"),
        "normalization_weight":       [],
        "normalization_weight_inclusive": [],
        "normalized_isr_weight":      get_shifts("isr"),
        "normalized_fsr_weight":      get_shifts("fsr"),
        # "normalized_pu_weight": get_shifts("minbias_xs"),
        # "normalized_njet_btag_deepjet_weight": get_shifts(*(f"btag_{u}" for u in cfg.x.btag_unc_names)),
        # "electron_weight": get_shifts("e"),
        # "muon_weight": get_shifts("mu"),
        # "tau_weight": get_shifts(*(f"tau_{u}" for u in cfg.x.tau_unc_names)),
        # "trigger_weight": get_shifts(*(f"trigger_{leg}" for leg in trigger_legs)),
    }
    # Store in the config (preserving DotDict interface)
    cfg.x.event_weights = DotDict(base_event_weights)
    
    # --- Per-dataset customizations ---
    for dataset in cfg.datasets:
        # Initialize empty mapping for each dataset (inherits from global)
        dataset.x.event_weights = {}
        if dataset.has_tag("ttbar"):
            dataset.x.event_weights["top_pt_weight"] = get_shifts("top_pt")
        elif dataset.has_tag("dy"):
            # Placeholder for Drell–Yan reweighting uncertainties
            dataset.x.event_weights["dy_weight"] = []  # TODO: add DY uncertainty sources
    
    # -----------------------------------------------------------------------------
    # Shift groups (for plotting / uncertainty grouping)
    # -----------------------------------------------------------------------------
    cfg.x.shift_groups = {
        "jec": _names_from_tag(("jec", "jer")),
        "tec": _names_from_tag("tec"),
        "eec": _names_from_tag(("ees", "eer")),
        "ees": _names_from_tag("ees"),
        "eer": _names_from_tag("eer"),
        "lepton_sf": [s.name for s in (*get_shifts("e"), *get_shifts("mu"))],
        "btag_sf": [s.name for s in get_shifts(*(f"btag_{u}" for u in cfg.x.btag_unc_names))],
        "pdf": [s.name for s in get_shifts("pdf")],
        "murmuf": [s.name for s in get_shifts("murmuf")],
        "pu": [s.name for s in get_shifts("minbias_xs")],
    }

    ################################################################################################
    # add external configs: channels, categories, met filters, triggers, variables
    ################################################################################################

    # channels
    cfg.add_channel(name="etau", id=1, label=r"$e\tau_{h}$")
    cfg.add_channel(name="mutau", id=2, label=r"$\mu\tau_{h}$")
    cfg.add_channel(name="tautau", id=3, label=r"$\tau_{h}\tau_{h}$")
    cfg.add_channel(name="ee", id=4, label=r"$ee$")
    cfg.add_channel(name="mumu", id=5, label=r"$\mu\mu$")
    cfg.add_channel(name="emu", id=6, label=r"$e\mu$")
    cfg.add_channel(name="c3e", id=14, label=r"$eee$")
    cfg.add_channel(name="c2emu", id=15, label=r"$ee\mu$")
    cfg.add_channel(name="ce2mu", id=16, label=r"$e\mu\mu$")
    cfg.add_channel(name="c3mu", id=17, label=r"$\mu\mu\mu$")
    cfg.add_channel(name="c4e", id=18, label=r"$eeee$")
    cfg.add_channel(name="c3emu", id=19, label=r"$eee\mu$")
    cfg.add_channel(name="c2e2mu", id=20, label=r"$ee\mu\mu$")
    cfg.add_channel(name="ce3mu", id=21, label=r"$e\mu\mu\mu$")
    cfg.add_channel(name="c4mu", id=22, label=r"$\mu\mu\mu\mu$")
    # to be implemented
    cfg.add_channel(name="c3etau", id=23, label=r"$eee\tau_{h}$")
    cfg.add_channel(name="c2emutau", id=24, label=r"$ee\mu\tau_{h}$")
    cfg.add_channel(name="ce2mutau", id=25, label=r"$e\mu\mu\tau_{h}$")
    cfg.add_channel(name="c3mutau", id=26, label=r"$\mu\mu\mu\tau_{h}$")
    cfg.add_channel(name="c2e2tau", id=27, label=r"$ee\tau_{h}\tau_{h}$")
    cfg.add_channel(name="cemu2tau", id=28, label=r"$e\mu\tau_{h}\tau_{h}$")
    cfg.add_channel(name="c2mu2tau", id=29, label=r"$\mu\mu\tau_{h}\tau_{h}$")
    cfg.add_channel(name="ce3tau", id=30, label=r"$e\tau_{h}\tau_{h}\tau_{h}$")
    cfg.add_channel(name="cmu3tau", id=31, label=r"$\mu\tau_{h}\tau_{h}\tau_{h}$")
    cfg.add_channel(name="c4tau", id=32, label=r"$\tau_{h}\tau_{h}\tau_{h}\tau_{h}$")
    cfg.add_channel(name="c2e0or1tau", id=33, label=r"$ee\  \leq 1\,\tau_{h}$")
    cfg.add_channel(name="cemu0or1tau", id=34, label=r"$e\mu\ \leq 1\,\tau_{h}$")
    cfg.add_channel(name="c2mu0or1tau", id=35, label=r"$\mu\mu\ \leq 1\,\tau_{h}$")

    # add categories
    from multilepton.config.categories import add_categories
    from multilepton.config.variables import add_variables
    from multilepton.config.met_filters import add_met_filters
    
    add_categories(cfg)
    add_variables(cfg)
    add_met_filters(cfg)
    add_triggers(cfg)
            
    return cfg
