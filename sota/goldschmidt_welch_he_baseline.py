"""Goldschmidt-division HE baseline for the arithmetic stage of Welch's test.

This module is deliberately *not* labelled as an implementation of a complete
MHE/MPC paper protocol. It is the reproducible HE-only arithmetic baseline used
to compare Goldschmidt division with the proposed inverse-square-root method.
"""

from dataclasses import dataclass
from typing import Tuple

from pp_test.engine.HEengine import HEEngine
from hedata.data import HEData
from operators.operator import HEOperator


@dataclass(frozen=True)
class PublicWelchBounds:
    """Public bounds required by the two Goldschmidt quotient evaluations.

    Values must be agreed before encryption. Deriving a lower bound from a
    decrypted dataset would invalidate a private evaluation; plaintext checks
    belong only in the experiment's post-run validation.
    """

    v_lower: float
    v_upper: float
    d_lower: float
    d_upper: float
    t2_upper: float
    df_upper: float

    def validate(self) -> None:
        for name, lower, upper in (
            ("V", self.v_lower, self.v_upper),
            ("D", self.d_lower, self.d_upper),
        ):
            if lower <= 0.0 or upper <= lower:
                raise ValueError(
                    f"{name} bounds must satisfy 0 < lower < upper; "
                    f"received ({lower}, {upper})."
                )
        if self.t2_upper <= 0.0 or self.df_upper <= 0.0:
            raise ValueError("t2_upper and df_upper must be positive.")


class GoldschmidtWelchHEBaseline:
    """HE-only Welch ``T^2, df`` baseline using normalized Goldschmidt division.

    It returns ``T^2`` rather than signed ``T`` because its two quotients are
    ``T^2=(mean1-mean2)^2/V`` and ``df=V^2/D``. The public bounds are explicit
    so the comparison remains auditable.
    """

    protocol_scope = "HE-only arithmetic-stage Goldschmidt baseline"

    def __init__(self, engine: HEEngine, verbose: bool = True):
        self._engine = engine
        self._ho = HEOperator(engine)
        self._verbose = verbose

    def bootstrap_count(self) -> int:
        return self._ho.bootstrap_count()

    def reset_bootstrap_count(self) -> None:
        self._ho.reset_bootstrap_count()

    def _goldschmidt_from_initial_factor(
        self,
        numerator: HEData,
        denominator: HEData,
        initial_factor: HEData,
        num_iter: int,
        bootstrap_level: int,
        bootstrap_threshold: int,
    ) -> HEData:
        """Run Goldschmidt in residual form with depth one per iteration.

        Let ``q_0 = a*w_0`` and ``r_0 = 1-b*w_0``.  The recurrence is

        ``q_(i+1) = q_i * (1 + r_i), r_(i+1) = r_i^2``.

        The two ciphertext multiplications in a round are independent, hence
        they have multiplicative depth one (although this Python reference
        invokes the evaluator calls sequentially).  Keeping only ``r`` avoids
        separately bootstrapping the algebraically redundant pair
        ``b_i`` and ``w_i = 2 - b_i`` used by the earlier prototype.
        """

        if num_iter < 1:
            raise ValueError("num_iter must be at least one.")

        ho = self._ho
        # q_0 and r_0 each consume one level and are independent.
        quotient = ho.mult(numerator, initial_factor)
        residual = ho.add_const(
            ho.mult_const(ho.mult(denominator, initial_factor), -1.0),
            1.0,
        )

        def refresh_state(keep_residual: bool) -> None:
            nonlocal quotient, residual
            if quotient.level() < bootstrap_threshold:
                if self._verbose:
                    print("[Goldschmidt] bootstrapping quotient")
                quotient = ho.do_bootstrapping(quotient, bootstrap_level)
            # The final residual is not consumed.  Bootstrapping it would add
            # work without improving the returned quotient.
            if keep_residual and residual.level() < bootstrap_threshold:
                if self._verbose:
                    print("[Goldschmidt] bootstrapping residual")
                residual = ho.do_bootstrapping(residual, bootstrap_level)

        if self._verbose:
            print(
                "[Goldschmidt] "
                f"iter=0, q_level={quotient.level()}, "
                f"r_level={residual.level()}"
            )
        refresh_state(keep_residual=num_iter > 1)

        for iteration in range(1, num_iter):
            # Both multiplications have the same depth increment of one.
            correction = ho.add_const(residual, 1.0)
            quotient = ho.mult(quotient, correction)
            residual = ho.mult(residual, residual)
            if self._verbose:
                print(
                    "[Goldschmidt] "
                    f"iter={iteration}, q_level={quotient.level()}, "
                    f"r_level={residual.level()}"
                )
            refresh_state(keep_residual=iteration < num_iter - 1)

        return quotient

    def HE_HDiv_Normalized(
        self,
        a: HEData,
        b: HEData,
        b_l: float,
        b_u: float,
        q_max: float,
        num_iter: int = 6,
        bootstrap_level: int = 11,
        bootstrap_threshold: int = 4,
        return_denormalized: bool = True,
    ) -> HEData:
        """Approximate ``a / b`` with normalized Goldschmidt iteration."""

        if b_l <= 0.0 or b_u <= b_l:
            raise ValueError("HDiv requires public bounds 0 < b_l < b_u.")
        if q_max <= 0.0:
            raise ValueError("q_max must be positive.")

        ho = self._ho
        b_norm = ho.mult_const(b, 1.0 / b_u)
        a_norm = ho.mult_const(a, 1.0 / (b_u * q_max))
        b_lower_norm = b_l / b_u

        # Linear initial approximation of 1 / b_norm on [b_lower_norm, 1].
        denominator = b_lower_norm**2 + 6.0 * b_lower_norm + 1.0
        w = ho.mult_const(b_norm, -8.0 / denominator)
        w = ho.add_const(w, 8.0 * (b_lower_norm + 1.0) / denominator)

        quotient = self._goldschmidt_from_initial_factor(
            a_norm,
            b_norm,
            w,
            num_iter,
            bootstrap_level,
            bootstrap_threshold,
        )

        if return_denormalized:
            return ho.mult_const(quotient, q_max)
        return quotient

    def HE_HDiv(
        self,
        a: HEData,
        b: HEData,
        b_l: float,
        b_u: float,
        num_iter: int = 6,
        bootstrap_level: int = 11,
        bootstrap_threshold: int = 4,
    ) -> HEData:
        """Legacy unnormalized Goldschmidt division API.

        New experiments should prefer :meth:`HE_HDiv_Normalized`, whose
        quotient bound makes the bootstrap range explicit.
        """

        if b_l <= 0.0 or b_u <= b_l:
            raise ValueError("HDiv requires public bounds 0 < b_l < b_u.")
        ho = self._ho
        denominator = b_l**2 + 6.0 * b_l * b_u + b_u**2
        w = ho.add_const(
            ho.mult_const(b, -8.0 / denominator),
            8.0 * (b_l + b_u) / denominator,
        )
        return self._goldschmidt_from_initial_factor(
            a,
            b,
            w,
            num_iter,
            bootstrap_level,
            bootstrap_threshold,
        )

    def welch_t2_df(
        self,
        x1: HEData,
        x2: HEData,
        bounds: PublicWelchBounds,
        num_iter: int = 6,
        return_debug: bool = False,
    ):
        """Compute encrypted Welch ``T^2`` and Welch--Satterthwaite ``df``."""

        bounds.validate()
        n1, n2 = x1.size(), x2.size()
        if n1 <= 1 or n2 <= 1:
            raise ValueError("Welch test requires both groups to have n > 1.")
        if self._verbose:
            print(f"n1, n2: {n1} {n2}")

        ho = self._ho
        s1, s2 = ho.sum(x1), ho.sum(x2)
        q1 = ho.sum(ho.mult(x1, x1))
        q2 = ho.sum(ho.mult(x2, x2))
        mean1 = ho.mult_const(s1, 1.0 / n1)
        mean2 = ho.mult_const(s2, 1.0 / n2)
        mean_difference = ho.sub(mean1, mean2)
        mean_difference_squared = ho.mult(mean_difference, mean_difference)
        variance1 = ho.mult_const(
            ho.sub(q1, ho.mult_const(ho.mult(s1, s1), 1.0 / n1)),
            1.0 / (n1 - 1),
        )
        variance2 = ho.mult_const(
            ho.sub(q2, ho.mult_const(ho.mult(s2, s2), 1.0 / n2)),
            1.0 / (n2 - 1),
        )
        a = ho.mult_const(variance1, 1.0 / n1)
        b = ho.mult_const(variance2, 1.0 / n2)
        v = ho.add(a, b)
        d = ho.add(
            ho.mult_const(ho.mult(a, a), 1.0 / (n1 - 1)),
            ho.mult_const(ho.mult(b, b), 1.0 / (n2 - 1)),
        )

        if self._verbose:
            print("[Goldschmidt] T2 = mean_difference_squared / V")
        t2 = self.HE_HDiv_Normalized(
            mean_difference_squared, v, bounds.v_lower, bounds.v_upper,
            bounds.t2_upper, num_iter=num_iter,
        )
        if self._verbose:
            print("[Goldschmidt] df = V_squared / D")
        df = self.HE_HDiv_Normalized(
            ho.mult(v, v), d, bounds.d_lower, bounds.d_upper,
            bounds.df_upper, num_iter=num_iter,
        )

        if return_debug:
            return {
                "T2": t2, "DF": df, "mean1": mean1, "mean2": mean2,
                "variance1": variance1, "variance2": variance2,
                "V": v, "D": d, "bounds": bounds,
                "protocol_scope": self.protocol_scope,
            }
        return t2, df

    # Backward-compatible name used by the original local prototype.
    def HE_Welch_T2_DF_HDiv_Baseline(
        self,
        x1: HEData,
        x2: HEData,
        V_bounds: Tuple[float, float],
        D_bounds: Tuple[float, float],
        T2_max: float = None,
        DF_max: float = None,
        R=None,
        num_iter: int = 6,
        return_debug: bool = False,
    ):
        if T2_max is None:
            if R is None:
                raise ValueError("T2_max or public input bound R is required.")
            T2_max = R**2 / V_bounds[0]
        if DF_max is None:
            DF_max = V_bounds[1] ** 2 / D_bounds[0]
        return self.welch_t2_df(
            x1, x2,
            PublicWelchBounds(
                V_bounds[0], V_bounds[1], D_bounds[0], D_bounds[1],
                T2_max, DF_max,
            ),
            num_iter=num_iter,
            return_debug=return_debug,
        )

    # Compatibility with the former secure_bioinformatics_baseline module.
    def HE_Welch_T2_DF(
        self,
        x1: HEData,
        x2: HEData,
        V_bounds: Tuple[float, float],
        D_bounds: Tuple[float, float],
        T2_max: float,
        DF_max: float,
        num_iter: int = 6,
        return_debug: bool = False,
    ):
        return self.HE_Welch_T2_DF_HDiv_Baseline(
            x1, x2, V_bounds, D_bounds, T2_max, DF_max,
            num_iter=num_iter,
            return_debug=return_debug,
        )
