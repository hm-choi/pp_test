from statistics import NormalDist

import numpy as np

from engine.HEengine import HEengine
from engine.HEdata import Ciphertext
from operators.HEapprox import HEApprox

from coeffs.t_critical_coeffs import (
    _T_CRITICAL_SQ_COEFFS,
    _T_CRITICAL_U_MIN,
    _T_CRITICAL_U_MAX,
)


# Lower bound of the normalized Welch variance V / V_max (public assumption),
# i.e. the inverse square root input domain is [V_BAR_MIN, 1].
V_BAR_MIN = 1e-2

# Inverse square root configuration.
#
# "select" (optional) maps the HE-DAP input level (level of the input before
# the domain mapping) to (log_degree, iteration, pre_bts), e.g. HE-DAP
# results; otherwise log_degree / iteration / pre_bts are used.
DEFAULT_INVSQRT = {
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

    @staticmethod
    def _check_public(
        n1: int,
        n2: int,
        R: float,
        test: str,
    ):

        if n1 <= 1 or n2 <= 1:
            raise ValueError(f"{test} requires both groups to have n > 1")

        if not np.isfinite(R) or R <= 0.0:
            raise ValueError("R must be a finite positive public bound")

    # ------------------------------------------------------------------
    # Public bounds and plaintext domain check
    # ------------------------------------------------------------------

    @staticmethod
    def welch_v_max(
        n1: int,
        n2: int,
        R: float,
    ):

        # For x in a public interval [L, U] with R = U - L, the unbiased
        # sample variance satisfies
        #
        # s_i^2 <= (n_i / (n_i - 1)) * R^2 / 4,
        #
        # so V = s1^2 / n1 + s2^2 / n2 <= R^2 / 4 * (1/(n1-1) + 1/(n2-1)).
        return (R ** 2 / 4.0) * (1.0 / (n1 - 1) + 1.0 / (n2 - 1))

    @staticmethod
    def variance_max(
        n: int,
        R: float,
    ):

        # Public upper bound of the unbiased sample variance for x in [L, U],
        # R = U - L.
        return n / (n - 1) * R ** 2 / 4.0

    @staticmethod
    def welch_d_max(
        n1: int,
        n2: int,
        R: float,
    ):

        # A_i = s_i^2 / n_i <= R^2 / (4 (n_i - 1)), so
        # D = A1^2 / (n1 - 1) + A2^2 / (n2 - 1)
        #   <= R^4 / 16 * (1/(n1-1)^3 + 1/(n2-1)^3).
        return (R ** 4 / 16.0) * (1.0 / (n1 - 1) ** 3 + 1.0 / (n2 - 1) ** 3)

    @classmethod
    def check_inv_sqrt_domain(
        cls,
        x1,
        x2,
        R: float,
    ):

        # Plaintext check of the inverse square root input domain of Welch's
        # t-test, V / V_max in [V_BAR_MIN, 1], run on the data before
        # encryption.
        #
        # Returns V / V_max; raises ValueError outside the domain.
        x1 = np.asarray(x1, dtype=np.float64)
        x2 = np.asarray(x2, dtype=np.float64)
        n1, n2 = len(x1), len(x2)

        value = (np.var(x1, ddof=1) / n1 + np.var(x2, ddof=1) / n2) / cls.welch_v_max(n1, n2, R)

        if not V_BAR_MIN <= value <= 1.0:
            raise ValueError(
                f"V / V_max = {value:.6g} outside [{V_BAR_MIN:g}, 1]"
            )

        return float(value)

    # ------------------------------------------------------------------
    # Shared building blocks
    # ------------------------------------------------------------------

    def _critical_value_from_inv_df(
        self,
        inv_df: Ciphertext,
        alpha: float,
        log_degree: int,
        min_output_level=0,
    ):

        # P(z) ~ t_{1-alpha/2, 1/u}^2 (squared two-sided t critical value) on
        # the normalized InvDF z. The input is bootstrapped beforehand when
        # needed so that the output keeps at least `min_output_level` levels.
        degree = 2 ** log_degree - 1
        table = _T_CRITICAL_SQ_COEFFS

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

        # InvDF = 1 / df with df >= 1, so InvDF lies in [0, 1].
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

        # The Chebyshev evaluation consumes log_degree levels.
        if z.level() < log_degree + min_output_level:
            self.approx._bootstrap(z)
        else:
            self.approx._ensure_cheb_level(z, degree)

        coeffs = self.engine._make_cheb_coeffs(coeffs_arr)

        return self.engine.evaluate_chebyshev(z, coeffs)

    @staticmethod
    def critical_value_sq_plain(
        u: float,
        alpha: float,
        log_degree: int,
    ):

        # Plaintext evaluation of the critical-value polynomial at u = 1 / df.
        coeffs = _T_CRITICAL_SQ_COEFFS[2 ** log_degree - 1][alpha]
        z = 2.0 * (u - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0

        return float(np.polynomial.chebyshev.chebval(z, coeffs))

    def _HE_sample_variance(
        self,
        x: Ciphertext,
        n: int,
    ):

        # mean = Sum(x) / n, Q = Sum(x * x) and the unbiased sample variance
        # s^2 = Q / (n - 1) - mean^2 * n / (n - 1).
        mean = self.engine.mult(self.engine.sum(x), 1.0 / n)
        q = self.engine.sum(self.engine.mult(x, x))

        var = self.engine.sub(
            self.engine.mult(q, 1.0 / (n - 1)),
            self.engine.mult(self.engine.mult(mean, mean), n / (n - 1)),
        )

        return mean, var

    def _inv_sqrt(
        self,
        v_bar: Ciphertext,
        min_level: int,
    ):

        # Approximate 1 / sqrt(V_bar) for V_bar in [V_BAR_MIN, 1]. The output
        # (|y| <= 1 / sqrt(V_BAR_MIN)) is refreshed with the extended
        # bootstrap when its level is below `min_level`.
        lo, hi = V_BAR_MIN, 1.0

        # First input for Newton iteration: V_bar / 2.
        x_half = self.engine.mult(v_bar, 0.5)

        # Second input for the Chebyshev approximation:
        # map V_bar from [lo, hi] to [-1, 1].
        x_cheb = self.engine.mult(v_bar, 2.0 / (hi - lo))
        x_cheb = self.engine.sub(x_cheb, (hi + lo) / (hi - lo))

        cfg = self.invsqrt

        if "select" in cfg:
            # Keyed by the HE-DAP input level l, i.e. the level of V_bar
            # before the one-level mapping (x_cheb is at level l - 1).
            log_degree, iteration, pre_bts = cfg["select"][x_cheb.level() + 1]
        else:
            log_degree, iteration, pre_bts = cfg["log_degree"], cfg["iteration"], cfg["pre_bts"]

        y = self.approx.invSqrt(
            x_half,
            x_cheb,
            domain=(lo, hi),
            log_degree=log_degree,
            iteration=iteration,
            pre_bts=pre_bts,
        )

        self.approx.ensure_level(y, min_level, 1.0 / np.sqrt(lo))

        return y

    def _step_from_score(
        self,
        normalized_score: Ciphertext,
        trace=None,
    ):

        # Approximate the sign of the normalized decision score (|score| < 1)
        # and map the sign output {-1, +1} to {0, 1}.
        sign_score = self.approx.sign(normalized_score)
        self._trace(trace, "sign", sign_score)

        step = self.engine.add(sign_score, 1.0)
        step = self.engine.mult(step, 0.5)

        return step

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

        # Encrypted quantities shared by every significance level:
        #
        # t      = T, the Welch t statistic,
        # inv_df = 1 / df (Welch-Satterthwaite),
        #
        # together with the public values used by the decision.
        self._check_public(n1, n2, R, "Welch t-test")

        v_max = self.welch_v_max(n1, n2, R)
        d_max = self.welch_d_max(n1, n2, R)

        # --------------------------------------------------------------
        # Welch t statistic
        # --------------------------------------------------------------

        # Means, unbiased sample variances and A_i = s_i^2 / n_i.
        mean1, var1 = self._HE_sample_variance(x1, n1)
        mean2, var2 = self._HE_sample_variance(x2, n2)

        a1 = self.engine.mult(var1, 1.0 / n1)
        a2 = self.engine.mult(var2, 1.0 / n2)

        # V = A1 + A2, normalized by the public bound: V_bar = V / V_max.
        v = self.engine.add(a1, a2)
        v_bar = self.engine.mult(v, 1.0 / v_max)

        # Y = 1 / sqrt(V_bar). Level 9 is needed for 1/V_bar, 1/V_bar^2,
        # 1/df and the critical-value mapping (bootstrapped at level >= 3
        # before the Chebyshev evaluation).
        y = self._inv_sqrt(v_bar, 9)
        self._trace(trace, "inv_sqrt_v_bar", y)

        # T = (mean1 - mean2) / sqrt(V_max) * Y.
        u = self.engine.mult(self.engine.sub(mean1, mean2), 1.0 / np.sqrt(v_max))
        t = self.engine.mult(u, y)
        self._trace(trace, "t", t)

        # --------------------------------------------------------------
        # Welch-Satterthwaite degrees of freedom
        # --------------------------------------------------------------

        # D = A1^2 / (n1 - 1) + A2^2 / (n2 - 1), normalized by the public
        # bound: D_bar = D / D_max in [0, 1]. The constants 1 / (n_i - 1) and
        # 1 / D_max are applied as one constant per term (one level).
        d1 = self.engine.mult(self.engine.mult(a1, a1), 1.0 / ((n1 - 1) * d_max))
        d2 = self.engine.mult(self.engine.mult(a2, a2), 1.0 / ((n2 - 1) * d_max))
        d_bar = self.engine.add(d1, d2)
        self.approx.ensure_level(d_bar, 7, 1.0)

        # 1 / V_bar = Y^2, 1 / V_bar^2 and
        # 1 / df = D / V^2 = D_bar / V_bar^2 * D_max / V_max^2.
        inv_v_bar = self.engine.mult(y, y)
        inv_v_bar_sq = self.engine.mult(inv_v_bar, inv_v_bar)

        inv_df = self.engine.mult(
            self.engine.mult(d_bar, inv_v_bar_sq),
            d_max / v_max ** 2,
        )
        self._trace(trace, "inv_df", inv_df)

        return {
            "t": t,
            "inv_df": inv_df,
            "R": R,
            "v_max": v_max,
            "n1": n1,
            "n2": n2,
        }

    def welch_score_bound(
        self,
        R: float,
        v_max: float,
        alpha: float,
        n1: int,
        n2: int,
        critical_log_degree=4,
    ):

        # Public bound of |s| with s = T^2 - c^2:
        # T^2 = (mean1 - mean2)^2 / V <= R^2 / (V_BAR_MIN * V_max) and
        # c^2 <= c_max^2. The Welch-Satterthwaite df is at least
        # min(n1, n2) - 1, so c_max^2 is the critical-value polynomial at
        # u = 1 / (min(n1, n2) - 1) (the polynomials of degree >= 15 are
        # positive and increasing in u on [0, 1]).
        c_max_sq = self.critical_value_sq_plain(
            1.0 / (min(n1, n2) - 1), alpha, critical_log_degree,
        )

        return max(R ** 2 / (V_BAR_MIN * v_max), c_max_sq)

    def welch_score(
        self,
        stats,
        alpha,
        critical_log_degree=4,
        trace=None,
    ):

        # Returns s / B (encrypted) with the decision score s = T^2 - c^2,
        # positive exactly when the null hypothesis is rejected, and the
        # public bound B = welch_score_bound(...).
        bound = self.welch_score_bound(
            stats["R"], stats["v_max"], alpha, stats["n1"], stats["n2"],
            critical_log_degree,
        )

        if trace is not None:
            trace["score_bound"] = bound

        # c^2 keeps 5 levels: T^2 - c^2 is not bootstrappable before the
        # normalization (|T^2 - c^2| <= B), and multiplying by a small 1 / B
        # takes up to two levels before StepHE (input level >= 3).
        critical_sq = self._critical_value_from_inv_df(
            stats["inv_df"],
            alpha,
            critical_log_degree,
            min_output_level=5,
        )
        self._trace(trace, "critical_sq", critical_sq)

        t_sq = self.engine.mult(stats["t"], stats["t"])
        self._trace(trace, "t_sq", t_sq)

        # s / B = (T^2 - c^2) / B.
        normalized_score = self.engine.mult(
            self.engine.sub(t_sq, critical_sq),
            1.0 / bound,
        )
        self._trace(trace, "score_normalized", normalized_score)

        return normalized_score

    def HE_Welch_decision(
        self,
        stats,
        alpha=0.05,
        critical_log_degree=4,
        trace=None,
    ):

        normalized_score = self.welch_score(
            stats, alpha, critical_log_degree, trace,
        )

        return self._step_from_score(normalized_score, trace)

    def HE_Welch_T_Test(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        R: float,
        alpha=0.05,
        critical_log_degree=4,
        trace=None,
    ):

        stats = self.HE_Welch_statistics(x1, x2, n1, n2, R, trace)

        return self.HE_Welch_decision(stats, alpha, critical_log_degree, trace)

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
    ):

        # Encrypted unbiased sample variances s_1^2, s_2^2 with the public
        # values used by the decision.
        self._check_public(n1, n2, R, "F-test")

        return {
            "var1": self._HE_sample_variance(x1, n1)[1],
            "var2": self._HE_sample_variance(x2, n2)[1],
            "n1": n1,
            "n2": n2,
            "R": R,
        }

    @classmethod
    def f_score_bound(
        cls,
        n1: int,
        n2: int,
        R: float,
        lower_critical: float,
        upper_critical: float,
    ):

        # Public bound of the decision score
        #
        # s = (s1^2 - F_L s2^2)(s1^2 - F_U s2^2):
        #
        # with s1^2 <= A, s2^2 <= B_v (public),
        # |s1^2 - F s2^2| <= max(A, F B_v), so |s| <= m_L * m_U.
        # Returns (m_L, m_U).
        var1_max = cls.variance_max(n1, R)
        var2_max = cls.variance_max(n2, R)

        return (
            max(var1_max, lower_critical * var2_max),
            max(var1_max, upper_critical * var2_max),
        )

    def f_score(
        self,
        stats,
        lower_critical: float,
        upper_critical: float,
        trace=None,
    ):

        # Returns s / B (encrypted) with the decision score
        #
        # s = (s1^2 - F_L s2^2)(s1^2 - F_U s2^2) = s2^4 (F - F_L)(F - F_U),
        #
        # which is positive exactly when F < F_L or F > F_U (s2^2 > 0), and
        # the public bound B = m_L * m_U (f_score_bound). Each factor is
        # normalized by its own bound.
        if (
            not np.isfinite(lower_critical)
            or not np.isfinite(upper_critical)
            or lower_critical <= 0.0
            or upper_critical <= lower_critical
        ):
            raise ValueError("Invalid public F critical values")

        n1, n2 = stats["n1"], stats["n2"]
        m_lower, m_upper = self.f_score_bound(
            n1, n2, stats["R"], lower_critical, upper_critical,
        )

        if trace is not None:
            trace["score_bound"] = m_lower * m_upper

        def factor(critical, m):
            # (s1^2 - critical * s2^2) / m.
            b = self.engine.sub(stats["var1"], self.engine.mult(stats["var2"], critical))
            return self.engine.mult(b, 1.0 / m)

        normalized_score = self.engine.mult(
            factor(lower_critical, m_lower),
            factor(upper_critical, m_upper),
        )
        self._trace(trace, "score_normalized", normalized_score)

        return normalized_score

    def HE_F_decision(
        self,
        stats,
        lower_critical: float,
        upper_critical: float,
        trace=None,
    ):

        # Two-sided F-test decision (reject when F < F_L or F > F_U).
        normalized_score = self.f_score(
            stats, lower_critical, upper_critical, trace,
        )

        return self._step_from_score(normalized_score, trace)

    def HE_F_Test(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        R: float,
        lower_critical: float,
        upper_critical: float,
        trace=None,
    ):

        stats = self.HE_F_statistics(x1, x2, n1, n2, R)

        return self.HE_F_decision(stats, lower_critical, upper_critical, trace)

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
        R: float,
        delta0=0.0,
        trace=None,
    ):

        # Encrypted Z = ((mean1 - mean2) - delta0) / sqrt(V_public) with the
        # public values used by the decision.
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

        if not np.isfinite(R) or R <= 0.0:
            raise ValueError("R must be a finite positive public bound")

        if not np.isfinite(delta0):
            raise ValueError("delta0 must be finite")

        mean1 = self.engine.mult(self.engine.sum(x1), 1.0 / n1)
        mean2 = self.engine.mult(self.engine.sum(x2), 1.0 / n2)

        numerator = self.engine.sub(mean1, mean2)

        if delta0 != 0.0:
            numerator = self.engine.sub(numerator, delta0)

        # The variance of the difference is public in the Z-test.
        v_public = sigma1_sq / n1 + sigma2_sq / n2

        z = self.engine.mult(numerator, 1.0 / np.sqrt(v_public))
        self._trace(trace, "z", z)

        return {
            "z": z,
            "v_public": v_public,
            # |(mean1 - mean2) - delta0| <= R + |delta0|.
            "diff_max": R + abs(delta0),
        }

    @staticmethod
    def z_score_bound(
        stats,
        alpha: float,
    ):

        # Public bound of |Z^2 - z^2|: Z^2 <= diff_max^2 / V_public, so
        # |Z^2 - z^2| <= max(diff_max^2 / V_public, z^2).
        critical_value = NormalDist().inv_cdf(1.0 - alpha / 2.0)

        return max(stats["diff_max"] ** 2 / stats["v_public"], critical_value ** 2)

    def HE_Z_decision(
        self,
        stats,
        alpha: float,
        trace=None,
    ):

        if not np.isfinite(alpha) or not 0.0 < alpha < 1.0:
            raise ValueError("alpha must be a finite value in (0, 1)")

        bound = self.z_score_bound(stats, alpha)

        if trace is not None:
            trace["score_bound"] = bound

        # Two-sided Gaussian critical value.
        critical_value = NormalDist().inv_cdf(1.0 - alpha / 2.0)

        # Compare Z^2 with the squared two-sided critical value.
        z_squared = self.engine.mult(stats["z"], stats["z"])
        score = self.engine.sub(z_squared, critical_value ** 2)
        self._trace(trace, "score", score)

        normalized_score = self.engine.mult(score, 1.0 / bound)

        return self._step_from_score(normalized_score, trace)

    def HE_Z_Test(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        sigma1_sq: float,
        sigma2_sq: float,
        R: float,
        delta0=0.0,
        alpha=None,
        trace=None,
    ):

        # Without alpha, returns the encrypted Z statistic.
        stats = self.HE_Z_statistic(
            x1, x2, n1, n2, sigma1_sq, sigma2_sq, R, delta0, trace,
        )

        if alpha is None:
            return stats["z"]

        return self.HE_Z_decision(stats, alpha, trace)
