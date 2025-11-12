import luigi

from columnflow.tasks.plotting import PlotVariables1D as CFPlotVariables1D

from multilepton.tasks.base import MultileptonTask


class PlotVariables1D(MultileptonTask, CFPlotVariables1D):
    """
    Wrapper around ColumnFlow's PlotVariables1D to include the MultileptonTask parameters.
    """
    task_namespace = "cf"
    
    def run(self):
        self.logger.info(f"Running PlotVariables1D with limit_dataset_files = {self.limit_dataset_files}")
        # Pass it to the config factory
        for config_name, factory in self.analysis.configs.items():
            factory_kwargs = {"limit_dataset_files": self.limit_dataset_files}
            config_obj = factory(**factory_kwargs)
        
        # Now pass the command-line parameter to the config factories
        for module, attr, name, cid in datasets:
            add_lazy_config(
                campaign_module=module,
                campaign_attr=attr,
                config_name=name,
                config_id=cid,
                add_limited=False,
                limit_dataset_files=self.limit_dataset_files
            )

        ## Continue with normal CF run
        super().run()
