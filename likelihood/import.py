"""Import models from MotifCentral, JASPAR, and HOCOMOCO."""
import json
import re

import numpy as np
import pyprobound
import requests
import torch


def import_motif_central(fit_id: int, checkpoint_name: str) -> None:
    model = json.loads(
        requests.get(
            "https://prod-gateway.motifcentral.org/"
            f"cellx/api/web/utility/fit/{fit_id}",
            timeout=5,
        ).text
    )

    alphabet = pyprobound.alphabets.DNA()
    kernel_size = int(model["modelSettings"]["bindingModes"][-1]["size"])
    pairwise_distance = min(
        kernel_size - 1,
        int(
            model["modelSettings"]["bindingModes"][-1]["dinucleotideDistance"]
        ),
    )
    psam = pyprobound.layers.PSAM(
        kernel_size=kernel_size,
        pairwise_distance=pairwise_distance,
        alphabet=alphabet,
        normalize=False,
    )

    mono = model["coefficients"]["bindingModes"][-1]["mononucleotide"]
    di = model["coefficients"]["bindingModes"][-1]["dinucleotide"]
    alphalen = len(alphabet.alphabet)
    for key, param in psam.betas.items():
        elements = [int(i) for i in key.split("-")]
        elements = [i - 1 for i in elements[:-1]] + elements[-1:]
        if len(elements) == 2:
            torch.nn.init.constant_(
                param, mono[elements[0] * alphalen + elements[1]]
            )
        else:
            torch.nn.init.constant_(
                param,
                di[elements[1] - elements[0] - 1][
                    elements[0] * (alphalen**2) + elements[2]
                ],
            )

    psam.save(checkpoint_name, [(0, 0)])


def import_hocomoco(species: str, fit_id: str, checkpoint_name: str) -> None:
    text = requests.get(
        "https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/"
        f"{species}/mono/pwm/{fit_id}.pwm",
        timeout=5,
    ).text
    array = [[float(j) for j in i.split("\t")] for i in text.split("\n")[1:]]
    array = [[j - max(i) for j in i] for i in array]

    alphabet = pyprobound.alphabets.DNA()
    psam = pyprobound.layers.PSAM(
        kernel_size=len(array), alphabet=alphabet, normalize=False
    )

    for key, param in psam.betas.items():
        elements = [int(i) for i in key.split("-")]
        elements = [i - 1 for i in elements[:-1]] + elements[-1:]
        if len(elements) == 2:
            torch.nn.init.constant_(param, array[elements[0]][elements[1]])

    psam.save(checkpoint_name, [(0, 0)])


def import_jaspar(fit_id: str, checkpoint_name: str) -> None:
    text = requests.get(
        f"https://jaspar.elixir.no/api/v1/matrix/{fit_id}.jaspar",
        timeout=5,
    ).text
    array = [
        [int(j) + 1 for j in re.split(r" |A|C|G|T|\[|\]", i) if len(j) > 0]
        for i in text.split("\n")[1:-1]
    ]
    array = list(map(list, zip(*array)))
    array = [[np.log(j) - np.log(max(i)) for j in i] for i in array]

    alphabet = pyprobound.alphabets.DNA()
    psam = pyprobound.layers.PSAM(
        kernel_size=len(array), alphabet=alphabet, normalize=False
    )

    for key, param in psam.betas.items():
        elements = [int(i) for i in key.split("-")]
        elements = [i - 1 for i in elements[:-1]] + elements[-1:]
        if len(elements) == 2:
            torch.nn.init.constant_(param, array[elements[0]][elements[1]])

    psam.save(checkpoint_name, [(0, 0)])


if __name__ == "__main__":
    for fit in [12715, 13165, 13774]:
        import_motif_central(
            fit, f"data/bindingModels/motifcentral_fit_{fit}.pt"
        )

    import_hocomoco(
        "HUMAN",
        "CTCF_HUMAN.H11MO.0.A",
        "data/bindingModels/HOCOMOCO_CTCF_HUMAN.H11MO.0.A.pt",
    )

    import_jaspar(
        "MA0139",
        "data/bindingModels/JASPAR_CTCF_MA0139.pt",
    )
