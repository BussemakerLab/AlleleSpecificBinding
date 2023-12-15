"""Allele-specific binding analysis."""
import dataclasses
from collections.abc import Iterable, Iterator
from typing import Literal, TypeVar, cast

import pandas as pd
import torch
from pandas import DataFrame
from pyprobound import (
    Aggregate,
    Contribution,
    CountBatch,
    Loss,
    LossModule,
    Table,
    __precision__,
)
from pyprobound.alphabets import Alphabet
from pyprobound.containers import TModuleList
from pyprobound.utils import betaln
from torch import Tensor
from typing_extensions import override

T = TypeVar("T", int, Tensor)


class ASBBatch(CountBatch):
    r"""A protocol for a set of rows from an ASB table.

    Attributes:
        seqs: A reference sequence tensor of shape
            :math:`(\text{minibatch},\text{length})` or
            :math:`(\text{minibatch},\text{in_channels},\text{length})`.
        alt_seqs: An alternate sequence tensor of shape
            :math:`(\text{minibatch},\text{length})` or
            :math:`(\text{minibatch},\text{in_channels},\text{length})`.
        target: A reference count tensor of shape
            :math:`(\text{minibatch},1)`.
        alt_count: An alternate count tensor of shape
            :math:`(\text{minibatch},1)`.
    """
    seqs: Tensor
    alt_seqs: Tensor
    target: Tensor
    alt_count: Tensor


@dataclasses.dataclass
class ASBBatchTuple(ASBBatch):
    r"""A NamedTuple for a set of rows from an ASB table.

    Attributes:
        seqs: A reference sequence tensor of shape
            :math:`(\text{minibatch},\text{length})` or
            :math:`(\text{minibatch},\text{in_channels},\text{length})`.
        alt_seqs: An alternate sequence tensor of shape
            :math:`(\text{minibatch},\text{length})` or
            :math:`(\text{minibatch},\text{in_channels},\text{length})`.
        target: A reference count tensor of shape
            :math:`(\text{minibatch},1)`.
        alt_count: An alternate count tensor of shape
            :math:`(\text{minibatch},1)`.
    """
    seqs: Tensor
    alt_seqs: Tensor
    target: Tensor
    alt_count: Tensor


def get_asb_dataframe(path: str) -> DataFrame:
    """Loads a tab-delimited ASB table into a Pandas dataframe."""
    return pd.read_csv(
        path,
        header=None,
        sep="\t",
        names=("left_flank", "ref", "alt", "right_flank", "ref_c", "alt_c"),
    )


class ASBTable(Table[ASBBatch]):
    """A tensor encoding of an ASB table with flank management.

    Attributes:
        left_flank_length (int): The length of the prepended sequence.
        right_flank_length (int): The length of the appended sequence.
    """

    def __init__(
        self,
        dataframe: DataFrame,
        alphabet: Alphabet,
        left_flank_length: int = 0,
        right_flank_length: int = 0,
        max_left_flank_length: int | None = None,
        max_right_flank_length: int | None = None,
    ) -> None:
        r"""Initializes the ASB table.

        Args:
            dataframe: The dataframe used to initialize the ASB table.
            alphabet: The alphabet used to encode sequences into tensors.
            left_flank_length: The initial length of the prepended sequence.
            right_flank_length: The initial length of the appended sequence.
            max_left_flank_length: The maximum allowed length of the prepended
                sequence.
            max_right_flank_length: The maximum allowed length of the appended
                sequence.
        """

        # Instance attributes
        self.alphabet = alphabet
        self.left_flank_length = 0
        self.right_flank_length = 0
        if max_left_flank_length is None:
            max_left_flank_length = dataframe["left_flank"].str.len().max()
        if max_right_flank_length is None:
            max_right_flank_length = dataframe["right_flank"].str.len().max()
        self.max_left_flank_length = max_left_flank_length
        self.max_right_flank_length = max_right_flank_length

        # Get dataframe data
        self.target = torch.tensor(
            dataframe["ref_c"].to_numpy(), dtype=__precision__
        )
        self.alt_count = torch.tensor(
            dataframe["alt_c"].to_numpy(), dtype=__precision__
        )
        self.left_flank = torch.stack(
            [alphabet.translate(seq) for seq in dataframe["left_flank"]]
        )
        self.right_flank = torch.stack(
            [alphabet.translate(seq) for seq in dataframe["right_flank"]]
        )
        self.ref = torch.stack(
            [alphabet.translate(seq) for seq in dataframe["ref"]]
        )
        self.alt = torch.stack(
            [alphabet.translate(seq) for seq in dataframe["alt"]]
        )
        self.seqs = self.ref.clone()
        self.alt_seqs = self.alt.clone()

        # Get number of probes per round
        self.counts_per_round = torch.stack(
            [self.target.sum(), self.alt_count.sum()]
        )

        # Set flank length
        super().__init__(
            left_flank_length=left_flank_length,
            right_flank_length=right_flank_length,
        )

    @override
    @property
    def input_shape(self) -> int:
        return self.seqs.shape[-1]

    @override
    @property
    def min_read_length(self) -> int:
        return self.max_read_length

    @override
    @property
    def max_read_length(self) -> int:
        return self.left_flank_length + self.right_flank_length + 1

    @override
    def set_flank_length(self, left: int = 0, right: int = 0) -> None:
        if left > self.max_left_flank_length:
            raise ValueError(
                f"left flank length of {left} exceeds"
                f" max_left_flank_length of {self.max_left_flank_length}"
            )
        if right > self.max_right_flank_length:
            raise ValueError(
                f"right flank length of {right} exceeds"
                f" max_right_flank_length of {self.max_right_flank_length}"
            )
        if left < 0 or right < 0:
            raise ValueError("Flank lengths must be nonnegative")
        self.left_flank_length = left
        self.right_flank_length = right
        tensors = [self.ref, self.right_flank[:, :right]]
        if left > 0:
            tensors.insert(0, self.left_flank[:, -left:])
        self.seqs = torch.cat(tensors, dim=-1)
        tensors[1] = self.alt
        self.alt_seqs = torch.cat(tensors, dim=-1)

    @override
    def __getitem__(self, idx: int) -> ASBBatchTuple:
        return ASBBatchTuple(
            self.seqs[idx],
            self.alt_seqs[idx],
            self.target[idx],
            self.alt_count[idx],
        )

    @override
    def __len__(self) -> int:
        return len(self.seqs)

    @override
    def get_setup_string(self) -> str:
        return "\n".join(
            [
                f"\t\tLeft Flank Length: {self.left_flank_length}",
                f"\t\tRight Flank Length: {self.right_flank_length}",
            ]
        )


class ASBAggregate(Aggregate):
    r"""PyProBound Aggregate with experiment-specific counts and dispersion.

    Attributes:
        contributions (TModuleList[Contribution]): The Contributions making up
            the aggregate.
    """

    unfreezable = Literal[Aggregate.unfreezable, "rho"]

    def __init__(
        self,
        contributions: Iterable[Contribution],
        train_rho: bool = True,
        log_rho: float = -1.0,
        name: str = "",
    ) -> None:
        super().__init__(
            contributions=contributions,
            train_concentration=False,
            target_concentration=1,
            name=name,
        )
        self.train_rho = train_rho
        self.log_rho = torch.nn.Parameter(
            torch.tensor(log_rho, dtype=__precision__), requires_grad=train_rho
        )

    @override
    def unfreeze(self, parameter: unfreezable = "all") -> None:
        """Unfreeze the desired parameter"""
        if parameter in ("rho", "all"):
            if self.train_rho:
                self.log_rho.requires_grad_()
        if parameter != "rho":
            super().unfreeze(parameter)


class ASBLoss(LossModule[ASBBatch]):
    """Multitask optimization of multiple ASB tables with Beta-Binomial loss.

    Attributes:
        aggregates (TModuleList[ASBAggregate]): The samples to be jointly
            optimized.
    """

    def __init__(
        self,
        aggregates: Iterable[ASBAggregate],
        lambda_l2: float = 1e-4,
        exponential_bound: float = 40,
    ) -> None:
        super().__init__(name="")

        # Store aggregate attribute
        self.aggregates: TModuleList[ASBAggregate] = TModuleList(aggregates)

        # Store loss attributes
        self.lambda_l2 = lambda_l2
        self.exponential_bound = exponential_bound

    @override
    def components(self) -> Iterator[ASBAggregate]:
        return iter(self.aggregates)

    @override
    def get_setup_string(self) -> str:
        out: list[str] = []
        out.extend(
            [
                "### Regularization:",
                f"\t L2 Lambda: {self.lambda_l2}",
                f"\t Exponential Bound: {self.exponential_bound}",
            ]
        )

        out.append("\n### Aggregates:")
        for agg in self.aggregates:
            out.append(f"\tAggregate: {str(agg)}")

        out.append("\n### Binding Components:")
        for binding_idx, (binding, optim) in enumerate(
            self.optim_procedure().items()
        ):
            out.extend([f"\t Mode {binding_idx}: {binding}"])
            out.append("\t\tFound In")
            for ancestors in optim.ancestry:
                out.append(f"\t\t\t{ancestors[-2]}")

        return "\n".join(out)

    def regularization(self, aggregate: ASBAggregate) -> Tensor:
        """Calculates parameter regularization.

        Args:
            experiment: The aggregate containing parameters to be regularized.

        Returns:
            The regularization value as a scalar tensor.
        """

        # Get flattened parameter vector
        param_vec = torch.cat(
            [
                p.flatten()
                for p in aggregate.parameters()
                if (torch.all(torch.isfinite(p)))
            ]
        )
        regularization = torch.tensor(0.0, device=param_vec.device)

        # Add L2 regularization
        if self.lambda_l2 > 0:
            regularization += self.lambda_l2 * param_vec.square().sum()

        # Add exponential barrier
        if self.exponential_bound != float("inf"):
            regularization += torch.sum(
                torch.exp(param_vec - self.exponential_bound)
                + torch.exp(-param_vec - self.exponential_bound)
            )

        return regularization

    @override
    def forward(self, batch: Iterable[ASBBatch]) -> Loss:
        """Calculates multitask weighted BetaBinomial NLL and regularization.

        Args:
            batch: Iterable of ASB tables to calculate the loss against.

        Returns:
            A NamedTuple with attributes `negloglik` (the BetaBinomial NLL) and
            `regularization`, both as scalar tensors.
        """

        negloglik = torch.tensor(0, dtype=__precision__)
        regularization = torch.tensor(0, dtype=__precision__)

        try:
            for agg, sample in zip(self.aggregates, batch, strict=True):
                device = agg.log_target_concentration.device

                x = sample.target.to(device)
                n = x + sample.alt_count.to(device)
                ref_score = agg(sample.seqs.to(device))
                alt_score = agg(sample.alt_seqs.to(device))

                prob = torch.exp(
                    ref_score - torch.logaddexp(ref_score, alt_score)
                )
                rho = torch.exp(agg.log_rho)

                a = prob * ((1 - rho) / rho)
                b = (1 - prob) * ((1 - rho) / rho)
                loglik = (
                    -torch.log(n + 1)
                    - betaln(n - x + 1, x + 1)
                    + betaln(x + a, n - x + b)
                    - betaln(a, b)
                )

                negloglik += -torch.sum(loglik) / len(x)
                regularization += self.regularization(agg)

        except ValueError as e:
            raise ValueError(
                "Length of aggregates and batches may not match?"
            ) from e

        return Loss(cast(Tensor, negloglik), cast(Tensor, regularization))
