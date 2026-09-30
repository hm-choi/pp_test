from __future__ import annotations

from functools import lru_cache
from typing import Sequence

import heaan as hn
import numpy as np


Monomial = tuple[int, ...]


def next_power_of_two_sqrt(x: int) -> int:
    if x <= 0:
        raise ValueError("x must be a positive integer")

    # k = ceil(log2(x) / 2)
    bitlen = x.bit_length() - 1  # floor(log2(x))
    k = (bitlen + 1) // 2        # ceil(bitlen / 2)

    if (1 << (2 * k)) < x:
        k += 1

    return 1 << k


def _validate_baby_step(baby_step: int) -> None:
    if baby_step <= 0 or baby_step & (baby_step - 1):
        raise ValueError(
            f"baby_step must be a positive power of two, got {baby_step}"
        )


def _make_blocks(
    coeffs: Sequence[float],
    baby_step: int,
) -> list[list[float]]:
    """
    Rewrite

        P(x) = sum_n a_n T_n(x)

    into

        P(x) = sum_k B_k(x) T_{kB}(x),

    where B = baby_step and each B_k contains only
    T_0, ..., T_{B-1}.

    Uses

        T_{m+r} = 2 T_m T_r - T_{m-r}.
    """
    degree = len(coeffs) - 1
    num_blocks = degree // baby_step + 1

    residual = list(coeffs)
    blocks = [
        [0.0] * baby_step
        for _ in range(num_blocks)
    ]

    for k in range(num_blocks - 1, 0, -1):
        base = k * baby_step

        for r in range(baby_step - 1, 0, -1):
            n = base + r

            if n > degree:
                continue

            value = residual[n]

            blocks[k][r] = 2.0 * value
            residual[base - r] -= value
            residual[n] = 0.0

        blocks[k][0] = residual[base]
        residual[base] = 0.0

    blocks[0][:min(baby_step, degree + 1)] = \
        residual[:baby_step]

    return blocks


def _make_giant_step_expander(baby_step: int):
    """
    Return expand(i), which expresses T_{iB} using products of

        T_B, T_{2B}, T_{4B}, ...

    via

        T_{m+n} = 2 T_m T_n - T_{|m-n|}.
    """

    @lru_cache(maxsize=None)
    def expand(i: int) -> dict[Monomial, int]:
        if i == 0:
            return {(): 1}

        # Powers of two are primitive giant-step terms.
        if i & (i - 1) == 0:
            return {(baby_step * i,): 1}

        high = 1 << (i.bit_length() - 1)
        low = i - high

        result: dict[Monomial, int] = {}

        for lhs, lhs_weight in expand(high).items():
            for rhs, rhs_weight in expand(low).items():
                monomial = tuple(sorted(lhs + rhs))

                result[monomial] = (
                    result.get(monomial, 0)
                    + 2 * lhs_weight * rhs_weight
                )

        for monomial, weight in expand(high - low).items():
            result[monomial] = (
                result.get(monomial, 0) - weight
            )

        return {
            monomial: weight
            for monomial, weight in result.items()
            if weight
        }

    return expand


def _channel_order(
    baby_step: int,
    num_blocks: int,
) -> list[Monomial]:
    """
    HEaaN BSGS channel order:

        T0
        TB
        T2B
        T2B*TB
        T4B
        T4B*TB
        T4B*T2B
        T4B*T2B*TB
        ...
    """
    if num_blocks <= 1:
        return [()]

    num_levels = (num_blocks - 1).bit_length()

    primitives = [
        baby_step << level
        for level in range(num_levels)
    ]

    order: list[Monomial] = [()]

    for level, high in enumerate(primitives):
        if level == 0:
            order.append((high,))
            continue

        for mask in range(1 << level):
            monomial = [high]

            for j in range(level):
                if mask & (1 << j):
                    monomial.append(primitives[j])

            order.append(tuple(sorted(monomial)))

    return order


def transform_coeffs_for_bsgs(
    coeffs: Sequence[float],
    baby_step: int,
) -> list[float]:
    """
    Convert ordinary Chebyshev coefficients

        [a_0, ..., a_d]

    to HEaaN's BSGS coefficient layout.

    The degree is inferred as len(coeffs) - 1 and may be arbitrary;
    it does not need to be 2^k - 1.
    """
    if len(coeffs) == 0:
        raise ValueError("coeffs must be non-empty")

    _validate_baby_step(baby_step)

    blocks = _make_blocks(coeffs, baby_step)
    expand = _make_giant_step_expander(baby_step)

    channel_order = _channel_order(
        baby_step,
        len(blocks),
    )

    # Allocate directly in final channel order.
    channel_index = {
        monomial: i
        for i, monomial in enumerate(channel_order)
    }

    transformed = [
        0.0
        for _ in range(len(channel_order) * baby_step)
    ]

    for block_index, block in enumerate(blocks):
        for monomial, weight in expand(block_index).items():
            offset = channel_index[monomial] * baby_step

            for j, coefficient in enumerate(block):
                transformed[offset + j] += weight * coefficient

    # The padded BSGS layout is always at least as long as the
    # original coefficient vector. Entries beyond this point are zero
    # for the partial final channel.
    return transformed[:len(coeffs)]


def make_hn_cheb(
    coeffs: Sequence[float],
    baby_step: int,
) -> hn.math.approx.ChebyshevCoefficients:
    transformed = np.asarray(
        transform_coeffs_for_bsgs(coeffs, baby_step),
        dtype=np.float64,
    )

    return hn.math.approx.ChebyshevCoefficients(
        transformed,
        baby_step,
    )