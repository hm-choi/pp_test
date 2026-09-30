"""Plaintext reference calculations for two-sample Welch's t-test."""

import numpy as np
from scipy import stats


def basic_stats(x: np.ndarray) -> tuple[int, float, float]:
    """
    Return the sample count, sum, and squared-sum.

    These statistics are sufficient to reconstruct the sample mean
    and unbiased sample variance.
    """

    values = np.asarray(x, dtype=np.float64)

    if values.ndim != 1 or values.size < 2:
        raise ValueError(
            "Welch's t-test requires a one-dimensional sample of size >= 2."
        )

    return (
        len(values),
        float(np.sum(values)),
        float(np.sum(values * values)),
    )


def welch_t_test_from_samples(
    x1: np.ndarray,
    x2: np.ndarray,
) -> dict[str, float | int]:
    """Compute Welch's t-test directly from two plaintext samples."""

    n1, sum1, square_sum1 = basic_stats(x1)
    n2, sum2, square_sum2 = basic_stats(x2)

    mean1, mean2 = sum1 / n1, sum2 / n2

    # Unbiased sample variance:
    # var = (Q - S^2 / n) / (n - 1)
    var1 = (
        square_sum1 - (sum1 * sum1) / n1
    ) / (n1 - 1)

    var2 = (
        square_sum2 - (sum2 * sum2) / n2
    ) / (n2 - 1)

    # Squared standard error:
    # SE^2 = var1 / n1 + var2 / n2
    se2 = var1 / n1 + var2 / n2

    if se2 <= 0.0:
        raise ValueError(
            "Welch's t-test is undefined when both sample variances are zero."
        )

    # Welch t-statistic:
    # t = (mean1 - mean2) / SE
    t_stat = (mean1 - mean2) / np.sqrt(se2)

    # Welch-Satterthwaite approximation for the degrees of freedom
    df_denominator = (
        ((var1 / n1) ** 2) / (n1 - 1)
        + ((var2 / n2) ** 2) / (n2 - 1)
    )
    df = (se2 ** 2) / df_denominator

    # Two-sided p-value
    p_value = 2.0 * stats.t.sf(abs(t_stat), df)

    return {
        "n1": n1,
        "n2": n2,
        "mean1": float(mean1),
        "mean2": float(mean2),
        "var1": float(var1),
        "var2": float(var2),
        "se": float(np.sqrt(se2)),
        "t_stat": float(t_stat),
        "df": float(df),
        "p_value": float(p_value),
    }


def welch_t_test_from_sufficient_statistics(
    x1: np.ndarray,
    x2: np.ndarray,
) -> dict[str, float | int]:
    """
    Compute Welch's t-test using sufficient statistics.

    This form mirrors the HE implementation, where the test is computed
    from aggregate statistics instead of directly using raw samples.
    """

    n1, sum1, square_sum1 = basic_stats(x1)
    n2, sum2, square_sum2 = basic_stats(x2)

    mean1, mean2 = sum1 / n1, sum2 / n2

    # Unbiased sample variance reconstructed from:
    # S = sum(x_i), Q = sum(x_i^2)
    #
    # var = (Q - S^2 / n) / (n - 1)
    var1 = (
        square_sum1 - (sum1 * sum1) / n1
    ) / (n1 - 1)

    var2 = (
        square_sum2 - (sum2 * sum2) / n2
    ) / (n2 - 1)

    # Variance contributions to the squared standard error
    a = var1 / n1
    b = var2 / n2
    se2 = a + b

    if se2 <= 0.0:
        raise ValueError(
            "Welch's t-test is undefined when both sample variances are zero."
        )

    # Welch t-statistic
    t_stat = (mean1 - mean2) / np.sqrt(se2)

    # Welch-Satterthwaite approximation:
    #
    #        (a + b)^2
    # df = -------------------------
    #      a^2/(n1-1) + b^2/(n2-1)
    df_denominator = (
        (a * a) / (n1 - 1)
        + (b * b) / (n2 - 1)
    )
    df = (se2 * se2) / df_denominator

    # Two-sided p-value
    p_value = 2.0 * stats.t.sf(abs(t_stat), df)

    return {
        "n1": n1,
        "n2": n2,
        "S1": sum1,
        "S2": sum2,
        "Q1": square_sum1,
        "Q2": square_sum2,
        "mean1": float(mean1),
        "mean2": float(mean2),
        "var1": float(var1),
        "var2": float(var2),
        "t_stat": float(t_stat),
        "df": float(df),
        "p_value": float(p_value),
    }