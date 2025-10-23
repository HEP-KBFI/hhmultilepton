# coding: utf-8

"""
Configuration of the HH → multi-leptons analysis.
"""

from __future__ import annotations

import importlib

import order as od

from columnflow.util import DotDict

from multilepton.hist_hooks.blinding import add_hooks as add_blinding_hooks
from multilepton.hist_hooks.qcd import add_hooks as add_qcd_hooks
from multilepton.hist_hooks.binning import add_hooks as add_binning_hooks
from multilepton.tasks.base import MultileptonTask
from multilepton.config.configs_multilepton import add_config



analysis_multilepton = od.Analysis(
    name="analysis_multilepton",
    id=1,
)

# analysis-global versions
# (empty since we use the lookup from the law.cfg instead)
analysis_multilepton.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
analysis_multilepton.x.bash_sandboxes = [
    "$CF_BASE/sandboxes/cf.sh",
    "$CF_BASE/sandboxes/venv_columnar.sh",
    "$MULTILEPTON_BASE/sandboxes/venv_columnar_tf.sh",
]

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
analysis_multilepton.x.cmssw_sandboxes = [
# "$CF_BASE/sandboxes/cmssw_default.sh",
]

#=======================================
# analysis-wide groups and defaults
#=======================================
# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
analysis_multilepton.x.config_groups = {}

# named function hooks that can modify store_parts of task outputs if needed
analysis_multilepton.x.store_parts_modifiers = {}

#=======================================
# hist hooks
#=======================================
analysis_multilepton.x.hist_hooks = DotDict()
add_blinding_hooks(analysis_multilepton)
add_qcd_hooks(analysis_multilepton)
add_binning_hooks(analysis_multilepton)


def add_lazy_config(
    *,
    campaign_module: str,
    campaign_attr: str,
    config_name: str,
    config_id: int,
    add_limited: bool = True,
    **kwargs,
):
    def create_factory(
        config_id: int,
        config_name_postfix: str = "",
        limit_dataset_files: int | None = None,
    ):
        def factory(configs: od.UniqueObjectIndex):
            # import the campaign
            mod = importlib.import_module(campaign_module)
            campaign = getattr(mod, campaign_attr)
            
            # use the task parameter if available
            #current_task = MultileptonTask.get_current_task()
            #if current_task and current_task.limit_dataset_files != -1:
            #    limit_dataset_files = current_task.limit_dataset_files
            return add_config(
                analysis_multilepton,
                campaign.copy(),
                config_name=config_name + config_name_postfix,
                config_id=config_id,
                limit_dataset_files=limit_dataset_files,
                **kwargs,
            )
        return factory

    analysis_multilepton.configs.add_lazy_factory(config_name, create_factory(config_id))
    if add_limited:
        analysis_multilepton.configs.add_lazy_factory(f"{config_name}_limited", 
                create_factory(config_id + 200, "_limited", 1))

#=======================================
# Priavte uhh NanoAOD datasets configs
#=======================================
# 2022, preEE, v14
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2022_preEE_nano_uhh_v14",
    campaign_attr="campaign_run3_2022_preEE_nano_uhh_v14",
    config_name="22preEE_v14_private",
    config_id=5014,
    add_limited=False,
)
# 2022, postEE, v14
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2022_postEE_nano_uhh_v14",
    campaign_attr="campaign_run3_2022_postEE_nano_uhh_v14",
    config_name="22postEE_v14_private",
    config_id=6014,
    add_limited=False,
)
# 2022, preEE, v14
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2022_preEE_nano_uhh_v14",
    campaign_attr="campaign_run3_2022_preEE_nano_uhh_v14",
    config_name="22preEE_v14_private",
    config_id=5114,
    add_limited=False,
)
# 2023, preBPix, v14
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2023_preBPix_nano_uhh_v14",
    campaign_attr="campaign_run3_2023_preBPix_nano_uhh_v14",
    config_name="23preBPix_v14_private",
    config_id=7014,
    add_limited=False,
)
# 2023, postBPix, v14
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2023_postBPix_nano_uhh_v14",
    campaign_attr="campaign_run3_2023_postBPix_nano_uhh_v14",
    config_name="23postBPix_v14_private",
    config_id=8014,
    add_limited=False,
)

#=======================================
# Central NanoAOD datasets configs
#=======================================
# 2022, preEE, v12
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2022_preEE_nano_v12",
    campaign_attr="campaign_run3_2022_preEE_nano_v12",
    config_name="22preEE_v12_central",
    config_id=5012,
    add_limited=False,
)
# 2022, postEE, v12
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2022_postEE_nano_v12",
    campaign_attr="campaign_run3_2022_postEE_nano_v12",
    config_name="22postEE_v12_central",
    config_id=6012,
    add_limited=False,
)
# 2023, preBPix, v12
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2023_preBPix_nano_v12",
    campaign_attr="campaign_run3_2023_preBPix_nano_v12",
    config_name="23preBPix_v12_central",
    config_id=7012,
    add_limited=False,
)
# 2023, postBPix, v12
add_lazy_config(
    campaign_module="cmsdb.campaigns.run3_2023_postBPix_nano_v12",
    campaign_attr="campaign_run3_2023_postBPix_nano_v12",
    config_name="23postBPix_v12_central",
    config_id=8012,
    add_limited=False,
)
