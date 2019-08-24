import os.path
import pandas as pd
import re
import math

def formatPath(path):
    """ Formats the path string so it must contain forward slashes and doesn't include them at the end

    Requires the path dir """

    path = path.rstrip("/") if path.endswith('/') else path
    path = path.replace("\\", "") if path.find("\\") else path

    return path

def load_matrix(matrix_file):
    """
    If using a raw PFM then convert it to PWM - psuedocounts ignored here
    Assumes only one frequency count is being used at a time
    """
    if os.path.isfile(matrix_file):
        pass
    else:
        Exception("File doesn't exist")
    formatted_pth = formatPath(matrix_file)
    fn = formatted_pth.split("/")[-1:]
    tmp_fn = formatted_pth.replace(fn[0], "tmp_file.jaspar")
    # Prep the file first
    with open(matrix_file) as f, open(tmp_fn, 'w') as out:
        for line in f:
            if ">" not in line:
                out.write(line)
    # Strip the pesky extra spaces with sed!
    os.system(f"sed -i 's/  \+/ /g' {tmp_fn}")
    pfm_df = pd.read_csv(tmp_fn, sep=" ", lineterminator='\n', header=None)
    pfm_df_transpose = pfm_df.transpose()
    col_len = 1/(len(pfm_df_transpose.index))
    pfm_df = pfm_df.astype(float)
    # Obtain PWM
    PWM_Df = pfm_df.applymap(lambda x: math.log2((x/col_len)))
    PWM_Df.rename(index={0:'A',1:'C', 2:'G', 3:'T'}, inplace=True)

    return PWM_Df
