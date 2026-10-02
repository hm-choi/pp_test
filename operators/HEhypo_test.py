from statistics import NormalDist

import numpy as np

from engine.HEengine import HEengine
from engine.HEdata import Ciphertext
from operators.HEapprox import HEApprox

from coeffs.t_critical_coeffs import (
    _T_CRITICAL_COEFFS,
    _T_CRITICAL_SQ_COEFFS,
    _T_CRITICAL_U_MIN,
    _T_CRITICAL_U_MAX,
)


# Inverse square root configuration.
#
# method="raw"        : x = V in [v_min, v_max(R, n)]; Chebyshev coefficients
#                       are generated from the public bounds.
# method="normalized" : x = V / v_max(R, n) in [v_min, 1]; fixed coefficients
#                       from coeffs/invSqrt_coeffs_normalized.py.
#
# "select" (optional) maps the HE-DAP input level (level of x before the
# mapping) to (log_degree, iteration, pre_bts), e.g. HE-DAP results;
# otherwise log_degree / iteration / pre_bts are used.
DEFAULT_INVSQRT = {
    "method": "raw",
    "v_min": 1e-3,
    "log_degree": 6,
    "iteration": 7,
    "pre_bts": False,
}


class HEHypothesisTesting:

    def __init__(
        self,
        engine: HEengine,
        invsqrt=None,
    ):

        self.engine = engine
        self.approx = HEApprox(engine)
        self.invsqrt = {**DEFAULT_INVSQRT, **(invsqrt or {})}
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

    @staticmethod
    def _trace(
        trace,
        name,
        ctxt: Ciphertext,
    ):

        # Keep a copy of an intermediate ciphertext for later decryption.
        if trace is not None:
            trace[name] = Ciphertext(ctxt)

    def _critical_value_from_inv_df(
        self,
        inv_df: Ciphertext,
        alpha: float,
        log_degree: int,
        target="t",
        scale=1.0,
    ):

        # target="t"  : P(u) ~ t_{1-alpha/2, 1/u}
        # target="t2" : P(u) ~ t_{1-alpha/2, 1/u}^2
        # The returned value is scale * P(u) (scale is folded into the
        # coefficients, so it costs no level).
        degree = 2 ** log_degree - 1

        if target == "t":
            table = _T_CRITICAL_COEFFS
        elif target == "t2":
            table = _T_CRITICAL_SQ_COEFFS
        else:
            raise ValueError(f"Unsupported critical-value target: {target}")

        if not np.isfinite(alpha) or not 0.0 < alpha < 1.0:
            raise ValueError("alpha must be a finite value in (0, 1)")

        if degree not in table:
            supported = ", ".join(str(key) for key in sorted(table))
            raise ValueError(
                f"Unsupported critical-value degree={degree}. Supported: {supported}"
            )

        if alpha not in table[degree]:
            supported = ", ".join(
                str(key) for key in sorted(table[degree])
            )
            raise ValueError(
                f"Unsupported alpha={alpha}. Supported: {supported}"
            )

        coeffs_arr = np.asarray(
            table[degree][alpha],
            dtype=np.float64,
        )

        if len(coeffs_arr) < 2 or not np.all(np.isfinite(coeffs_arr)):
            raise ValueError("Critical-value coefficients must be finite")

        # InvDF = 1 / df with df in [1, 2000], so
        # InvDF lies in [1/2000, 1].
        u_min = _T_CRITICAL_U_MIN
        u_max = _T_CRITICAL_U_MAX

        # Map InvDF from [u_min, u_max] to the Chebyshev domain [-1, 1]:
        #
        # z = 2 * InvDF / (u_max - u_min)
        #     - (u_max + u_min) / (u_max - u_min).
        z = self.engine.mult(
            inv_df,
            2.0 / (u_max - u_min),
        )
        z = self.engine.sub(
            z,
            (u_max + u_min) / (u_max - u_min),
        )

        self.approx._ensure_cheb_level(z, degree)

        coeffs = self.engine._make_cheb_coeffs(coeffs_arr * scale)

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

    def _inv_sqrt_from_terms(
        self,
        terms,
        v_max: float,
        min_level: int,
    ):

        # Approximate 1 / sqrt(x) for
        #
        # x = x_scale * V,  V = sum_i factor_i * ctxt_i,
        #
        # using the configured domain method, and return (y, x_scale) with
        # y ~ 1 / sqrt(x). The factors are applied directly to ctxt_i to
        # avoid additional ciphertext-scalar multiplications after forming x;
        # callers fold x_scale into their own public constants instead of
        # spending a level on rescaling y.
        #
        # y is refreshed (extended bootstrap, |y| <= 1 / sqrt(lo)) when its
        # level is below `min_level`.
        cfg = self.invsqrt
        v_min = cfg["v_min"]

        if cfg["method"] == "normalized":
            # x = V / v_max in [v_min, 1].
            lo, hi = v_min, 1.0
            norm = 1.0 / v_max
        elif cfg["method"] == "raw":
            # x = V in [v_min, v_max].
            lo, hi = v_min, v_max
            norm = 1.0
        else:
            raise ValueError(f"Unsupported invSqrt method: {cfg['method']}")

        if hi <= lo:
            raise ValueError("v_max must be greater than v_min")

        def combine(scale, shift=0.0):
            ret = None

            for ctxt, factor in terms:
                part = self.engine.mult(ctxt, factor * norm * scale)
                ret = part if ret is None else self.engine.add(ret, part)

            if shift != 0.0:
                ret = self.engine.sub(ret, shift)

            return ret

        # First input for Newton iteration: x / 2.
        x_half = combine(0.5)

        # Second input for the Chebyshev approximation:
        # map x from [lo, hi] to [-1, 1].
        x_cheb = combine(2.0 / (hi - lo), (hi + lo) / (hi - lo))

        if "select" in cfg:
            # Keyed by the HE-DAP input level l, i.e. the level of x before
            # the one-level mapping (x_cheb is at level l - 1).
            log_degree, iteration, pre_bts = cfg["select"][x_cheb.level() + 1]
        else:
            log_degree, iteration, pre_bts = cfg["log_degree"], cfg["iteration"], cfg["pre_bts"]

        y = self.approx.invSqrt(
            x_half,
            x_cheb,
            log_degree=log_degree,
            iteration=iteration,
            domain=(lo, hi),
            method=cfg["method"],
            pre_bts=pre_bts,
        )

        self.approx.ensure_level(y, min_level, 1.0 / np.sqrt(lo))

        return y, norm

    # ------------------------------------------------------------------
    # Welch's t-test
    # ------------------------------------------------------------------

    def HE_Welch_statistics(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        R: float,
        trace=None,
    ):

        # Returns encrypted (T^2, 1/df); both are shared by every
        # significance level.
        if n1 <= 1 or n2 <= 1:
            raise ValueError("Welch t-test requires both groups to have n > 1")

        if not np.isfinite(R) or R <= 0.0:
            raise ValueError("R must be a finite positive public bound")

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
        v_max = (R ** 2 / 4.0) * (
            1.0 / (n1 - 1) + 1.0 / (n2 - 1)
        )

        # Compute encrypted sufficient statistics for the two variances.
        #
        # scaled_var_i = n_i * (n_i - 1) * s_i^2.
        scaled_var1 = self._HE_var_scaled(x1, n1)
        scaled_var2 = self._HE_var_scaled(x2, n2)

        # --------------------------------------------------------------
        # Inverse square root of V
        # --------------------------------------------------------------

        # a_i = s_i^2 / n_i = scaled_var_i / (n_i^2 * (n_i - 1)).
        #
        # inv_sqrt_x ~ 1 / sqrt(x) with x = x_scale * V. Level 7 is needed
        # for 1/x, 1/x^2, 1/df and the critical-value mapping (bootstrapped
        # at level >= 3 before the Chebyshev evaluation).
        inv_sqrt_x, x_scale = self._inv_sqrt_from_terms(
            [
                (scaled_var1, 1.0 / n1 / n1 / (n1 - 1)),
                (scaled_var2, 1.0 / n2 / n2 / (n2 - 1)),
            ],
            v_max,
            7,
        )
        self._trace(trace, "inv_sqrt_x", inv_sqrt_x)

        if trace is not None:
            trace["x_scale"] = x_scale

        # Square once to obtain 1 / x = 1 / (x_scale * V).
        inv_x = self.engine.mult(inv_sqrt_x, inv_sqrt_x)

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
        #
        # x_scale is also folded in (D * x_scale^2 / x^2 = D / V^2).
        scale1_for_df = x_scale / (n1 ** 2 * (n1 - 1) ** 1.5)
        scale2_for_df = x_scale / (n2 ** 2 * (n2 - 1) ** 1.5)

        d1_base = self.engine.mult(scaled_var1, scale1_for_df)
        d2_base = self.engine.mult(scaled_var2, scale2_for_df)

        # Squaring gives a_i^2 / (n_i - 1).
        d1 = self.engine.mult(d1_base, d1_base)
        d2 = self.engine.mult(d2_base, d2_base)

        d = self.engine.add(d1, d2)

        # Compute 1 / x^2. (No bootstrap here: 1 / x^2 is not bounded by 1;
        # the critical-value input z in [-1, 1] is bootstrapped instead.)
        inv_x_squared = self.engine.mult(inv_x, inv_x)

        # Welch-Satterthwaite inverse degree of freedom:
        #
        # 1 / df = D / V^2 = (D * x_scale^2) / x^2.
        inv_df = self.engine.mult(d, inv_x_squared)
        self._trace(trace, "inv_df", inv_df)

        # --------------------------------------------------------------
        # Welch t statistic
        # --------------------------------------------------------------

        # Compute the encrypted difference of sample means, scaled by
        # sqrt(x_scale).
        mean_scale = np.sqrt(x_scale)

        mean1 = self.engine.mult(
            self.engine.sum(x1),
            mean_scale / n1,
        )
        mean2 = self.engine.mult(
            self.engine.sum(x2),
            mean_scale / n2,
        )

        t_numerator = self.engine.sub(mean1, mean2)
        self._trace(trace, "mean_diff_scaled", t_numerator)

        t_numerator_squared = self.engine.mult(
            t_numerator,
            t_numerator,
        )

        # T^2 = (mean1 - mean2)^2 / V = (x_scale * (mean1 - mean2)^2) / x.
        t_squared = self.engine.mult(
            t_numerator_squared,
            inv_x,
        )
        self._trace(trace, "t_squared", t_squared)

        return t_squared, inv_df

    def HE_Welch_decision(
        self,
        t_squared: Ciphertext,
        inv_df: Ciphertext,
        alpha=0.05,
        critical_log_degree=4,
        score_bound=1.0,
        target="t2",
        trace=None,
    ):

        if not np.isfinite(score_bound) or score_bound <= 0.0:
            raise ValueError("score_bound must be a finite positive public bound")

        # Approximate the two-sided t critical value (or its square)
        # from 1 / df. The score normalization 1 / score_bound is folded
        # into the critical-value coefficients (1 / sqrt(score_bound) for
        # target "t", which is squared afterwards), because the Chebyshev
        # output level is too low for further scalar multiplications
        # before the sign approximation.
        #
        # trace: "critical_scaled" = critical / sqrt(B) (target "t") or
        #        critical^2 / B (target "t2"); "score_normalized" = score / B.
        if target == "t":
            critical_scale = 1.0 / np.sqrt(score_bound)
        else:
            critical_scale = 1.0 / score_bound

        critical_scaled = self._critical_value_from_inv_df(
            inv_df,
            alpha,
            critical_log_degree,
            target,
            critical_scale,
        )
        self._trace(trace, "critical_scaled", critical_scaled)

        if target == "t":
            critical_squared_scaled = self.engine.mult(
                critical_scaled,
                critical_scaled,
            )
        else:
            critical_squared_scaled = critical_scaled

        # A positive score means that T^2 exceeds the squared
        # two-sided critical value:
        #
        # score / score_bound = T^2 / score_bound - critical^2 / score_bound.
        t_squared_scaled = self.engine.mult(
            t_squared,
            1.0 / score_bound,
        )
        normalized_score = self.engine.sub(
            t_squared_scaled,
            critical_squared_scaled,
        )
        self._trace(trace, "score_normalized", normalized_score)

        # Approximate the sign of the decision score.
        sign_score = self.approx.sign(normalized_score)
        self._trace(trace, "sign", sign_score)

        # Map the sign output {-1, +1} to {0, 1}.
        step = self.engine.add(sign_score, 1.0)
        step = self.engine.mult(step, 0.5)

        return step

    def HE_Welch_T_Test(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        R: float,
        alpha=0.05,
        critical_log_degree=4,
        score_bound=1.0,
        target="t2",
        trace=None,
    ):

        t_squared, inv_df = self.HE_Welch_statistics(
            x1, x2, n1, n2, R, trace,
        )

        return self.HE_Welch_decision(
            t_squared,
            inv_df,
            alpha,
            critical_log_degree,
            score_bound,
            target,
            trace,
        )

    # ------------------------------------------------------------------
    # F-test
    # ------------------------------------------------------------------

    def HE_F_statistics(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        R: float,
        trace=None,
    ):

        # Returns encrypted (x_scale * variance1, 1 / (x_scale * variance2));
        # their product is F. Shared by every significance level.
        if n1 <= 1 or n2 <= 1:
            raise ValueError("F-test requires n1 > 1 and n2 > 1")

        if not np.isfinite(R) or R <= 0.0:
            raise ValueError("R must be a finite positive public bound")

        # For x in [0, R], the unbiased sample variance satisfies
        #
        # s^2 <= n / (n - 1) * R^2 / 4.
        var2_max = (
            n2
            / (n2 - 1)
            * R ** 2
            / 4.0
        )

        # scaled_var_i = n_i * (n_i - 1) * s_i^2.
        scaled_var1 = self._HE_var_scaled(x1, n1)
        scaled_var2 = self._HE_var_scaled(x2, n2)

        # --------------------------------------------------------------
        # Inverse square root of variance2
        # --------------------------------------------------------------

        # inv_sqrt_x ~ 1 / sqrt(x) with x = x_scale * variance2. Level 7 is
        # needed for 1/x, F, the decision score and its normalization
        # (sign bootstraps at level >= 3).
        inv_sqrt_x, x_scale = self._inv_sqrt_from_terms(
            [(scaled_var2, 1.0 / n2 / (n2 - 1))],
            var2_max,
            7,
        )
        self._trace(trace, "inv_sqrt_x", inv_sqrt_x)

        if trace is not None:
            trace["x_scale"] = x_scale

        # 1 / x = 1 / (x_scale * variance2).
        inv_x = self.engine.mult(
            inv_sqrt_x,
            inv_sqrt_x,
        )

        # variance1 * x_scale, so that (variance1 * x_scale) / x = F.
        variance1_scaled = self.engine.mult(
            scaled_var1,
            x_scale / n1 / (n1 - 1),
        )
        self._trace(trace, "variance1_scaled", variance1_scaled)

        return variance1_scaled, inv_x

    def HE_F_decision(
        self,
        variance1_scaled: Ciphertext,
        inv_x: Ciphertext,
        lower_critical: float,
        upper_critical: float,
        score_bound: float,
        trace=None,
    ):

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

        # F = variance1 / variance2.
        f_for_score = self.engine.mult(
            variance1_scaled,
            inv_x,
        )
        self._trace(trace, "f_statistic", f_for_score)

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
        self._trace(trace, "score", score)

        normalized_score = self.engine.mult(
            score,
            1.0 / score_bound,
        )

        sign_score = self.approx.sign(normalized_score)
        self._trace(trace, "sign", sign_score)

        # Map {-1, +1} to {0, 1}.
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
        trace=None,
    ):

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

        variance1_scaled, inv_x = self.HE_F_statistics(
            x1, x2, n1, n2, R, trace,
        )

        df1 = n1 - 1
        df2 = n2 - 1

        if lower_critical is None:
            # F = variance1 / variance2.
            f_statistic = self.engine.mult(
                variance1_scaled,
                inv_x,
            )

            return f_statistic, df1, df2

        return self.HE_F_decision(
            variance1_scaled,
            inv_x,
            lower_critical,
            upper_critical,
            score_bound,
            trace,
        )

    # ------------------------------------------------------------------
    # Z-test
    # ------------------------------------------------------------------

    def HE_Z_statistic(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        sigma1_sq: float,
        sigma2_sq: float,
        delta0=0.0,
        trace=None,
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
        self._trace(trace, "z", z)

        return z

    def HE_Z_decision(
        self,
        z: Ciphertext,
        alpha: float,
        score_bound: float,
        trace=None,
    ):

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
        self._trace(trace, "score", score)

        normalized_score = self.engine.mult(
            score,
            1.0 / score_bound,
        )

        sign_score = self.approx.sign(normalized_score)
        self._trace(trace, "sign", sign_score)

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
        trace=None,
    ):

        if alpha is None and score_bound is not None:
            raise ValueError(
                "score_bound requires an alpha for decision mode"
            )

        z = self.HE_Z_statistic(
            x1, x2, n1, n2, sigma1_sq, sigma2_sq, delta0, trace,
        )

        if alpha is None:
            return z

        return self.HE_Z_decision(z, alpha, score_bound, trace)
