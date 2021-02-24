import os
import pandas as pd

module_root_dir = os.path.split(__file__)[0]

results_specifications = pd.read_csv(os.path.join(module_root_dir, "resources", "results_specifications.csv"), sep='\t')
