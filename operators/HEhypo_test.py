from statistics import NormalDist

import numpy as np

from engine.HEengine import HEengine
from engine.HEdata import Ciphertext
from pp_test.operators.HEapprox import HEApprox

from coeffs.t_critical_coeffs import _T_CRITICAL_COEFFS


class HEHypothesisTesting:

    def __init__(
        self,
        engine: HEengine,
    ):

        self.engine = engine
        self.approx = HEApprox(engine)
        self._bootstrap_count = 0

    def bootstrap_count(
        self,
    ):

        return self._bootstrap_count + self.approx.bootstrap_count()

    def reset_bootstrap_count(
        self,
    ):

        self._bootstrap_count = 0
        self.approx.reset_bootstrap_count()

    def _bootstrap(
        self,
        ctxt: Ciphertext,
    ):

        # Bootstrap is performed in place.
        self.engine.bootstrap(ctxt)
        self._bootstrap_count += 1

    def _critical_value_from_inv_df(
        self,
        inv_df: Ciphertext,
        alpha: float,
        log_degree: int,
    ):

        degree = 2 ** log_degree - 1

        if not np.isfinite(alpha) or not 0.0 < alpha < 1.0:
            raise ValueError("alpha must be a finite value in (0, 1)")

        if degree not in _T_CRITICAL_COEFFS:
            supported = ", ".join(str(key) for key in sorted(_T_CRITICAL_COEFFS))
            raise ValueError(
                f"Unsupported critical-value degree={degree}. Supported: {supported}"
            )

        if alpha not in _T_CRITICAL_COEFFS[degree]:
            supported = ", ".join(
                str(key) for key in sorted(_T_CRITICAL_COEFFS[degree])
            )
            raise ValueError(
                f"Unsupported alpha={alpha}. Supported: {supported}"
            )

        coeffs_arr = np.asarray(
            _T_CRITICAL_COEFFS[degree][alpha],
            dtype=np.float64,
        )

        if len(coeffs_arr) < 2 or not np.all(np.isfinite(coeffs_arr)):
            raise ValueError("Critical-value coefficients must be finite")

        # Map InvDF from [0, 1] to the Chebyshev domain [-1, 1]:
        # z = 2 * InvDF - 1 = 2 / df - 1.
        z = self.engine.add(inv_df, inv_df)
        z = self.engine.sub(z, 1.0)

        # Ensure sufficient levels before Chebyshev evaluation.
        self.engine._ensure_cheb_level(z, degree)

        coeffs = self.engine._make_cheb_coeffs(coeffs_arr)

        return self.engine.evaluate_chebyshev(z, coeffs)

    def _HE_var_scaled(
        self,
        x: Ciphertext,
        n: int,
    ):

        # Compute
        #
        # n * sum(x_i^2) - (sum(x_i))^2
        #
        # which is equal to
        #
        # n * (n - 1) * sample_variance.
        s = self.engine.sum(x)
        q = self.engine.sum(self.engine.mult(x, x))

        nq = self.engine.mult(q, n)
        s_squared = self.engine.mult(s, s)

        return self.engine.sub(nq, s_squared)

    def HE_Welch_T_Test(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        R: float,
        alpha=0.05,
        critical_log_degree=15,
        score_bound=1.0,
    ):

        if n1 <= 1 or n2 <= 1:
            raise ValueError("Welch t-test requires both groups to have n > 1")

        if not np.isfinite(R) or R <= 0.0:
            raise ValueError("R must be a finite positive public bound")

        if not np.isfinite(score_bound) or score_bound <= 0.0:
            raise ValueError("score_bound must be a finite positive public bound")

        # Welch variance term:
        #
        # V = s1^2 / n1 + s2^2 / n2.
        #
        # For x in [0, R], the unbiased sample variance satisfies
        #
        # s_i^2 <= (n_i / (n_i - 1)) * R^2 / 4.
        #
        # Therefore,
        #
        # s_i^2 / n_i <= R^2 / (4 * (n_i - 1)),
        #
        # which gives the following public upper bound for V.
        # We assume V lies in [v_min, v_max].
        v_min = 1e-3
        v_max = (R ** 2 / 4.0) * (
            1.0 / (n1 - 1) + 1.0 / (n2 - 1)
        )

        if v_max <= v_min:
            raise ValueError("v_max must be greater than v_min")

        # Compute encrypted sufficient statistics for the two variances.
        #
        # scaled_var_i = n_i * (n_i - 1) * s_i^2.
        scaled_var1 = self._HE_var_scaled(x1, n1)
        scaled_var2 = self._HE_var_scaled(x2, n2)

        # --------------------------------------------------------------
        # Inverse square root of V
        # --------------------------------------------------------------

        # First input for Newton iteration:
        #
        # V / 2 = (a1 + a2) / 2,
        #
        # where
        #
        # a_i = s_i^2 / n_i.
        #
        # The factor is applied directly to scaled_var_i to avoid an
        # additional ciphertext-scalar multiplication after forming V.
        scale1_half = 1.0 / n1 / n1 / (n1 - 1) / 2.0
        scale2_half = 1.0 / n2 / n2 / (n2 - 1) / 2.0

        a1_half = self.engine.mult(scaled_var1, scale1_half)
        a2_half = self.engine.mult(scaled_var2, scale2_half)

        v_half = self.engine.add(a1_half, a2_half)

        # Second input for the Chebyshev approximation:
        #
        # map V from [v_min, v_max] to [-1, 1].
        #
        # V_cheb
        # = 2V / (v_max - v_min)
        #   - (v_max + v_min) / (v_max - v_min).
        #
        # Again, the multiplicative part of the mapping is applied
        # directly to scaled_var_i to reduce level consumption.
        scale1_cheb = (
            1.0
            / n1
            / n1
            / (n1 - 1)
            * 2.0
            / (v_max - v_min)
        )
        scale2_cheb = (
            1.0
            / n2
            / n2
            / (n2 - 1)
            * 2.0
            / (v_max - v_min)
        )

        a1_cheb = self.engine.mult(scaled_var1, scale1_cheb)
        a2_cheb = self.engine.mult(scaled_var2, scale2_cheb)

        v_cheb = self.engine.add(a1_cheb, a2_cheb)
        v_cheb = self.engine.sub(
            v_cheb,
            (v_max + v_min) / (v_max - v_min),
        )

        # Approximate 1 / sqrt(V).
        inv_sqrt_v = self.approx.invSqrt(v_half, v_cheb)

        # Square once to obtain 1 / V.
        inv_v = self.engine.mult(inv_sqrt_v, inv_sqrt_v)

        # --------------------------------------------------------------
        # Welch-Satterthwaite degrees of freedom
        # --------------------------------------------------------------

        # Define
        #
        # D = a1^2 / (n1 - 1) + a2^2 / (n2 - 1),
        #
        # where a_i = s_i^2 / n_i.
        #
        # Since
        #
        # scaled_var_i = n_i * (n_i - 1) * s_i^2,
        #
        # multiplying scaled_var_i by
        #
        # 1 / (n_i^2 * (n_i - 1)^(3/2))
        #
        # directly gives
        #
        # a_i / sqrt(n_i - 1).
        scale1_for_df = 1.0 / (n1 ** 2 * (n1 - 1) ** 1.5)
        scale2_for_df = 1.0 / (n2 ** 2 * (n2 - 1) ** 1.5)

        d1_base = self.engine.mult(scaled_var1, scale1_for_df)
        d2_base = self.engine.mult(scaled_var2, scale2_for_df)

        # Squaring gives a_i^2 / (n_i - 1).
        d1 = self.engine.mult(d1_base, d1_base)
        d2 = self.engine.mult(d2_base, d2_base)

        d = self.engine.add(d1, d2)

        # Compute 1 / V^2.
        inv_v_squared = self.engine.mult(inv_v, inv_v)
        self._bootstrap(inv_v_squared)

        # Welch-Satterthwaite inverse degree of freedom:
        #
        # 1 / df = D / V^2.
        inv_df = self.engine.mult(d, inv_v_squared)

        # Approximate the two-sided t critical value from 1 / df.
        critical_value = self._critical_value_from_inv_df(
            inv_df,
            alpha,
            critical_log_degree,
        )

        # --------------------------------------------------------------
        # Welch t statistic
        # --------------------------------------------------------------

        # Compute the encrypted difference of sample means.
        mean1 = self.engine.mult(
            self.engine.sum(x1),
            1.0 / n1,
        )
        mean2 = self.engine.mult(
            self.engine.sum(x2),
            1.0 / n2,
        )

        t_numerator = self.engine.sub(mean1, mean2)
        t_numerator_squared = self.engine.mult(
            t_numerator,
            t_numerator,
        )

        # Refresh a separate copy because inv_v is also used in the
        # degree-of-freedom branch above.
        inv_v_for_score = Ciphertext(inv_v)
        self._bootstrap(inv_v_for_score)

        # T^2 = (mean1 - mean2)^2 / V.
        t_squared = self.engine.mult(
            t_numerator_squared,
            inv_v_for_score,
        )

        # --------------------------------------------------------------
        # Hypothesis-test decision
        # --------------------------------------------------------------

        critical_squared = self.engine.mult(
            critical_value,
            critical_value,
        )

        # A positive score means that T^2 exceeds the squared
        # two-sided critical value.
        score = self.engine.sub(
            t_squared,
            critical_squared,
        )
        normalized_score = self.engine.mult(
            score,
            1.0 / score_bound,
        )

        # Approximate the sign of the decision score.
        sign_score = self.approx.sign(normalized_score)

        # Map the sign output {-1, +1} to {0, 1}.
        step = self.engine.add(sign_score, 1.0)
        step = self.engine.mult(step, 0.5)

        return step

    def HE_F_Test(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        R: float,
        lower_critical=None,
        upper_critical=None,
        score_bound=None,
    ):

        if n1 <= 1 or n2 <= 1:
            raise ValueError("F-test requires n1 > 1 and n2 > 1")

        if not np.isfinite(R) or R <= 0.0:
            raise ValueError("R must be a finite positive public bound")

        decision_args = (
            lower_critical,
            upper_critical,
            score_bound,
        )

        if any(value is not None for value in decision_args):
            if any(value is None for value in decision_args):
                raise ValueError(
                    "F decision mode requires lower_critical, upper_critical, "
                    "and score_bound together"
                )

            if (
                not np.isfinite(lower_critical)
                or not np.isfinite(upper_critical)
                or not np.isfinite(score_bound)
                or lower_critical <= 0.0
                or upper_critical <= lower_critical
                or score_bound <= 0.0
            ):
                raise ValueError(
                    "Invalid public F critical values or score_bound"
                )

        # For x in [0, R], the unbiased sample variance satisfies
        #
        # s^2 <= n / (n - 1) * R^2 / 4.
        #
        # A positive lower bound is required by the inverse square root
        # approximation used for the denominator variance.
        var_min = 1e-3
        var2_max = (
            n2
            / (n2 - 1)
            * R ** 2
            / 4.0
        )

        if var2_max <= var_min:
            raise ValueError("var2_max must be greater than var_min")

        # scaled_var_i = n_i * (n_i - 1) * s_i^2.
        scaled_var1 = self._HE_var_scaled(x1, n1)
        scaled_var2 = self._HE_var_scaled(x2, n2)

        # Recover variance1 directly.
        variance1 = self.engine.mult(
            scaled_var1,
            1.0 / n1 / (n1 - 1),
        )

        # --------------------------------------------------------------
        # Inverse square root of variance2
        # --------------------------------------------------------------

        # First input:
        #
        # variance2 / 2.
        scale2_half = (
            1.0
            / n2
            / (n2 - 1)
            / 2.0
        )

        var2_half = self.engine.mult(
            scaled_var2,
            scale2_half,
        )

        # Second input:
        #
        # map variance2 from [var_min, var2_max] to [-1, 1].
        scale2_cheb = (
            1.0
            / n2
            / (n2 - 1)
            * 2.0
            / (var2_max - var_min)
        )

        var2_cheb = self.engine.mult(
            scaled_var2,
            scale2_cheb,
        )
        var2_cheb = self.engine.sub(
            var2_cheb,
            (var2_max + var_min) / (var2_max - var_min),
        )

        # Approximate 1 / sqrt(variance2).
        inv_sqrt_var2 = self.approx.invSqrt(
            var2_half,
            var2_cheb,
        )

        # 1 / variance2.
        inv_var2 = self.engine.mult(
            inv_sqrt_var2,
            inv_sqrt_var2,
        )

        # F = variance1 / variance2.
        f_statistic = self.engine.mult(
            variance1,
            inv_var2,
        )

        df1 = n1 - 1
        df2 = n2 - 1

        if lower_critical is None:
            return f_statistic, df1, df2

        # Refresh a copy of the inverse denominator variance for
        # the decision path.
        inv_var2_for_score = Ciphertext(inv_var2)
        self._bootstrap(inv_var2_for_score)

        f_for_score = self.engine.mult(
            variance1,
            inv_var2_for_score,
        )

        # The product is positive when F lies outside
        # [lower_critical, upper_critical].
        lower_score = self.engine.sub(
            f_for_score,
            lower_critical,
        )
        upper_score = self.engine.sub(
            f_for_score,
            upper_critical,
        )
        score = self.engine.mult(
            lower_score,
            upper_score,
        )

        normalized_score = self.engine.mult(
            score,
            1.0 / score_bound,
        )

        sign_score = self.approx.sign(normalized_score)

        # Map {-1, +1} to {0, 1}.
        step = self.engine.add(sign_score, 1.0)
        step = self.engine.mult(step, 0.5)

        return step

    def HE_Z_Test(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        sigma1_sq: float,
        sigma2_sq: float,
        delta0=0.0,
        alpha=None,
        score_bound=None,
    ):

        if n1 <= 0 or n2 <= 0:
            raise ValueError("Z-test requires n1 > 0 and n2 > 0")

        if (
            not np.isfinite(sigma1_sq)
            or not np.isfinite(sigma2_sq)
            or sigma1_sq <= 0.0
            or sigma2_sq <= 0.0
        ):
            raise ValueError(
                "sigma1_sq and sigma2_sq must be finite positive values"
            )

        if not np.isfinite(delta0):
            raise ValueError("delta0 must be finite")

        if alpha is None and score_bound is not None:
            raise ValueError(
                "score_bound requires an alpha for decision mode"
            )

        if alpha is not None:
            if not np.isfinite(alpha) or not 0.0 < alpha < 1.0:
                raise ValueError(
                    "alpha must be a finite value in (0, 1)"
                )

            if (
                score_bound is None
                or not np.isfinite(score_bound)
                or score_bound <= 0.0
            ):
                raise ValueError(
                    "decision mode requires a finite positive score_bound"
                )

        # Compute the encrypted difference of sample means.
        mean1 = self.engine.mult(
            self.engine.sum(x1),
            1.0 / n1,
        )
        mean2 = self.engine.mult(
            self.engine.sum(x2),
            1.0 / n2,
        )

        numerator = self.engine.sub(
            mean1,
            mean2,
        )

        if delta0 != 0.0:
            numerator = self.engine.sub(
                numerator,
                delta0,
            )

        # The variance of the difference is public in the Z-test.
        v_public = (
            sigma1_sq / n1
            + sigma2_sq / n2
        )

        if v_public <= 0.0:
            raise ValueError("V_public must be positive")

        # Z = ((mean1 - mean2) - delta0) / sqrt(V_public).
        z = self.engine.mult(
            numerator,
            1.0 / np.sqrt(v_public),
        )

        if alpha is None:
            return z

        # Two-sided Gaussian critical value.
        critical_value = NormalDist().inv_cdf(
            1.0 - alpha / 2.0
        )

        # Compare Z^2 with the squared two-sided critical value.
        z_squared = self.engine.mult(z, z)
        score = self.engine.sub(
            z_squared,
            critical_value ** 2,
        )

        normalized_score = self.engine.mult(
            score,
            1.0 / score_bound,
        )

        sign_score = self.approx.sign(normalized_score)

        # Map {-1, +1} to {0, 1}.
        step = self.engine.add(sign_score, 1.0)
        step = self.engine.mult(step, 0.5)

        return step