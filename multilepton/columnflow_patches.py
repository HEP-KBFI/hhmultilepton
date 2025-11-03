# coding: utf-8

"""
Collection of patches of underlying columnflow tasks.
"""

import os
import law
import getpass

from columnflow.util import memoize


logger = law.logger.get_logger(__name__)


@memoize
def patch_columnar_pyarrow_version():
    """
    Comments out the pyarrow==21.0.0 line in the columnar.txt sandbox file.
    """
    columnar_path = os.path.join(
        os.environ["MULTILEPTON_BASE"], "modules", "columnflow", "sandboxes", "columnar.txt"
    )

    if not os.path.exists(columnar_path):
        logger.warning(f"File not found: {columnar_path}")
        return
    with open(columnar_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if "pyarrow==" in line and not line.strip().startswith("#"):
            new_lines.append(f"# {line.strip()}\n")
            logger.debug("Commented out pyarrow line in columnar.txt")
        else:
            new_lines.append(line)

    with open(columnar_path, "w") as f:
        f.writelines(new_lines)
    logger.info(f"Patched {columnar_path}: commented out pyarrow requirement")


@memoize
def patch_bundle_repo_exclude_files():
    """
    Patches the exclude_files attribute of the existing BundleRepo task to exclude files specific to _this_ analysis
    project.
    """
    from columnflow.tasks.framework.remote import BundleRepo

    # get the relative path to CF_BASE
    cf_rel = os.path.relpath(os.environ["CF_BASE"], os.environ["MULTILEPTON_BASE"])
    # amend exclude files to start with the relative path to CF_BASE
    exclude_files = [os.path.join(cf_rel, path) for path in BundleRepo.exclude_files]
    # add additional files
    exclude_files.extend([
        "docs", "tests", "data", "assets", ".law", ".setups", ".data", ".github",
    ])
    # overwrite them
    BundleRepo.exclude_files[:] = exclude_files
    logger.debug(f"patched exclude_files of {BundleRepo.task_family}")


@memoize
def patch_remote_workflow_poll_interval():
    """
    Patches the HTCondorWorkflow and SlurmWorkflow tasks to change the
    default value of the poll_interval parameter to 1 minute.
    """
    from columnflow.tasks.framework.remote import HTCondorWorkflow, SlurmWorkflow

    HTCondorWorkflow.poll_interval._default = 1.0  # minutes
    SlurmWorkflow.poll_interval._default = 1.0  # minutes
    logger.debug(f"patched poll_interval._default of {HTCondorWorkflow.task_family} and {SlurmWorkflow.task_family}")


@memoize
def patch_merge_reduction_stats_inputs():
    """
    Patches the MergeReductionStats task to set the default value of n_inputs to -1, so as to use all files to infer
    merging factors with full statistical precision.
    """
    from columnflow.tasks.reduction import MergeReductionStats

    MergeReductionStats.n_inputs._default = -1
    logger.debug(f"patched n_inputs default value of {MergeReductionStats.task_family}")


@memoize
def patch_htcondor_workflow_naf_resources():
    """
    Patches the HTCondorWorkflow task to declare user-specific resources when running on the NAF.
    """
    from columnflow.tasks.framework.remote import HTCondorWorkflow

    def htcondor_job_resources(self, job_num, branches):
        # one "naf_<username>" resource per job, indendent of the number of branches in the job
        return {f"naf_{getpass.getuser()}": 1}

    HTCondorWorkflow.htcondor_job_resources = htcondor_job_resources
    logger.debug(f"patched htcondor_job_resources of {HTCondorWorkflow.task_family}")


@memoize
def patch_slurm_partition_setting():
    """
    Patches the slurm remote workflow to allow setting things like partition
    by commandline instead of overiding with central default.
    """
    from columnflow.tasks.framework.remote import RemoteWorkflow
    
    RemoteWorkflow.exclude_params_branch.remove("slurm_partition")
    RemoteWorkflow.slurm_partition.significant = True
    RemoteWorkflow.exclude_params_branch.remove("slurm_flavor")
    RemoteWorkflow.slurm_flavor._choices.add("manivald")
    logger.debug(f"patched slurm partition/flavor settings of {RemoteWorkflow.task_family}")


@memoize
def patch_all():
    patch_bundle_repo_exclude_files()
    patch_remote_workflow_poll_interval()
    patch_slurm_partition_setting()
    patch_merge_reduction_stats_inputs()
    patch_columnar_pyarrow_version()    
    #patch_htcondor_workflow_naf_resources()
