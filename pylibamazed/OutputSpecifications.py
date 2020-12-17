import os
import pandas as pd

module_root_dir = os.path.split(__file__)[0]

candidate_specifications = pd.read_csv(os.path.join(module_root_dir, "resources", "candidates_specifications.csv"), sep='\t')
