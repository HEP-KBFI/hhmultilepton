# coding: utf-8

import law

from multilepton.columnflow_patches import patch_all

law.contrib.load("pandas")
# apply cf patches once
patch_all()
