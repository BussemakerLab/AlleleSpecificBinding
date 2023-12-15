"""Create ASB table for ASB-derived modeling."""
import pandas as pd

asb_df = pd.read_csv("../likelihood/data/accB/acc_undup_CTCF.csv")
ref_df = pd.read_csv("../likelihood/data/seq/CTCF_k30_ref.txt", header=None)
alt_df = pd.read_csv("../likelihood/data/seq/CTCF_k30_alt.txt", header=None)
new_df = pd.DataFrame(
    {
        "left": ref_df[0].str[:29],
        "ref": ref_df[0].str[29],
        "alt": alt_df[0].str[29],
        "right": ref_df[0].str[30:],
        "ref_c": asb_df["ref_c"],
        "alt_c": asb_df["alt_c"],
    }
)
new_df.to_csv("CTCF.tsv", sep="\t", header=False, index=False)
