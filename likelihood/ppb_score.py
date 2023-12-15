"""Score a PyProBound model from a checkpoint"""

import argparse

import pandas as pd
import pyprobound
import torch


def score_seqs(model: str, in_file: str, out_file: str):
    # Read input data
    df = pd.read_csv(in_file, header=None, index_col=0, sep="\t")
    df[1] = torch.ones((len(df),))
    df = df[[1]]

    # Convert to tensor dataset
    alphabet = pyprobound.alphabets.DNA()
    count_table = pyprobound.CountTable(df, alphabet)

    # Load model
    psam = pyprobound.layers.PSAM(kernel_size=1, alphabet=alphabet)
    conv1d = pyprobound.layers.Conv1d.from_psam(psam, count_table)
    mode = pyprobound.Mode([conv1d])
    psam.reload(model)

    with torch.inference_mode():
        df[1] = mode(count_table.seqs)

    df.to_csv(out_file, sep="\t", header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Score sequences using a PyProBound model"
    )
    parser.add_argument(
        "model",
        type=str,
        help="Path to PyProBound checkpoint",
    )
    parser.add_argument(
        "in_file",
        type=str,
        help="Path to read sequences from",
    )
    parser.add_argument(
        "out_file",
        type=str,
        help="Path to write output to",
    )

    args = parser.parse_args()
    score_seqs(model=args.model, in_file=args.in_file, out_file=args.out_file)
