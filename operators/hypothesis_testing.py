import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np
import heaan as hn

from engine.HEengine import HEengine
from engine.HEdata import Ciphertext
from operators.invSqrt import HEStats


class HEHypothesisTesting:

    _COEFFICIENT_DIR = Path(__file__).resolve().parents[1] / "coefficients"
    _T_CRITICAL_COEFFICIENTS = _COEFFICIENT_DIR / "t_dist_coeffs_by_degree.json"

    def __init__(self, engine: HEengine):

        self.engine = engine
        self.stats = HEStats(engine)
        self._bootstrap_count = 0

    def bootstrap_count(self):

        return self._bootstrap_count + self.stats.bootstrap_count()

    def reset_bootstrap_count(self):

        self._bootstrap_count = 0
        self.stats.reset_bootstrap_count()

    def _bootstrap(self, ctxt: Ciphertext):

        self.engine.bootstrap(ctxt)
        self._bootstrap_count += 1

    def _critical_value_from_inv_df(
        self,
        inv_df: Ciphertext,
        alpha: float,
        degree: int,
    ):

        if not np.isfinite(alpha) or not 0.0 < alpha < 1.0:
            raise ValueError("alpha must be a finite value in (0, 1)")

        path = self._T_CRITICAL_COEFFICIENTS

        if not path.is_file():
            raise FileNotFoundError(f"Missing t critical-value coefficients: {path}")

        with path.open("r", encoding="utf-8") as file:
            table = json.load(file)

        degree_key = f"degree_{degree}"
        alpha_key = f"alpha_{alpha:g}"

        if degree_key not in table:
            supported = ", ".join(key.removeprefix("degree_") for key in table)
            raise ValueError(
                f"Unsupported critical-value degree={degree}. Supported: {supported}"
            )

        if alpha_key not in table[degree_key]:
            supported = ", ".join(
                key.removeprefix("alpha_") for key in table[degree_key]
            )
            raise ValueError(f"Unsupported alpha={alpha}. Supported: {supported}")

        try:
            coefficients = np.asarray(
                [float(pair[0]) for pair in table[degree_key][alpha_key]],
                dtype=np.float64,
            )
        except (IndexError, TypeError, ValueError) as exc:
            raise ValueError("Invalid t critical-value coefficient format") from exc

        if len(coefficients) < 2 or not np.all(np.isfinite(coefficients)):
            raise ValueError("Critical-value coefficients must be finite")

        # z = 2 / df - 1 maps InvDF in [0, 1] to [-1, 1].
        z = self.engine.sub(self.engine.add(inv_df, inv_df), 1.0)

        required_levels = math.ceil(math.log2(len(coefficients) - 1)) + 3

        if z.level() < required_levels:
            self._bootstrap(z)

        cheb_coeffs = hn.math.approx.ChebyshevCoefficients(
            coefficients,
            len(coefficients),
        )

        return self.engine.evaluate_chebyshev(z, cheb_coeffs)

    def HE_Welch_T_Test(
        self,
        x1: Ciphertext,
        x2: Ciphertext,
        n1: int,
        n2: int,
        R: float,
        alpha=0.05,
        critical_degree=15,
        score_bound=1.0,
    ):

        if n1 <= 1 or n2 <= 1:
            raise ValueError("Welch t-test requires both groups to have n > 1")

        if not np.isfinite(R) or R <= 0.0:
            raise ValueError("R must be a finite positive public bound")

        if not np.isfinite(score_bound) or score_bound <= 0.0:
            raise ValueError("score_bound must be a finite positive public bound")

        s1 = self.engine.sum(x1)
        s2 = self.engine.sum(x2)

        q1 = self.engine.sum(self.engine.mult(x1, x1))
        q2 = self.engine.sum(self.engine.mult(x2, x2))

        mean1 = self.engine.mult(s1, 1.0 / n1)
        mean2 = self.engine.mult(s2, 1.0 / n2)

        s1_squared = self.engine.mult(s1, s1)
        s2_squared = self.engine.mult(s2, s2)

        variance1 = self.engine.sub(q1, self.engine.mult(s1_squared, 1.0 / n1))
        variance2 = self.engine.sub(q2, self.engine.mult(s2_squared, 1.0 / n2))

        variance1 = self.engine.mult(variance1, 1.0 / (n1 - 1))
        variance2 = self.engine.mult(variance2, 1.0 / (n2 - 1))

        a1 = self.engine.mult(variance1, 1.0 / n1)
        a2 = self.engine.mult(variance2, 1.0 / n2)
        v = self.engine.add(a1, a2)

        # For x in [0, R], sample variance is bounded by approximately R^2 / 4.
        v_max = (R ** 2 / 4.0) * (1.0 / (n1 - 1) + 1.0 / (n2 - 1))

        if v_max <= 0.0:
            raise ValueError("Computed V_max must be positive")

        v_norm = self.engine.mult(v, 1.0 / v_max)
        self._bootstrap(v_norm)

        inv_sqrt_v = self.stats.invSqrt(v_norm)

        t_numerator = self.engine.sub(mean1, mean2)
        t_numerator = self.engine.mult(t_numerator, 1.0 / np.sqrt(v_max))

        a1_squared = self.engine.mult(a1, a1)
        a2_squared = self.engine.mult(a2, a2)

        d1 = self.engine.mult(a1_squared, 1.0 / (n1 - 1))
        d2 = self.engine.mult(a2_squared, 1.0 / (n2 - 1))
        d = self.engine.add(d1, d2)

        a1_max = (R ** 2 / 4.0) / (n1 - 1)
        a2_max = (R ** 2 / 4.0) / (n2 - 1)
        d_max = a1_max ** 2 / (n1 - 1) + a2_max ** 2 / (n2 - 1)

        if d_max <= 0.0:
            raise ValueError("Computed D_max must be positive")

        d_norm = self.engine.mult(d, 1.0 / d_max)
        self._bootstrap(d_norm)

        # InvDF = D / V^2.
        inv_v_norm = self.engine.mult(inv_sqrt_v, inv_sqrt_v)
        inv_v_norm_squared = self.engine.mult(inv_v_norm, inv_v_norm)
        self._bootstrap(inv_v_norm_squared)

        inv_df = self.engine.mult(d_norm, inv_v_norm_squared)
        inv_df = self.engine.mult(inv_df, d_max / v_max ** 2)

        critical_value = self._critical_value_from_inv_df(
            inv_df,
            alpha,
            critical_degree,
        )

        t_numerator_squared = self.engine.mult(t_numerator, t_numerator)

        inv_v_for_score = Ciphertext(inv_v_norm)
        self._bootstrap(inv_v_for_score)

        t_squared = self.engine.mult(t_numerator_squared, inv_v_for_score)
        critical_squared = self.engine.mult(critical_value, critical_value)

        score = self.engine.sub(t_squared, critical_squared)
        normalized_score = self.engine.mult(score, 1.0 / score_bound)

        self._bootstrap(normalized_score)

        sign_score = self.stats.sign(normalized_score)

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

        decision_args = (lower_critical, upper_critical, score_bound)

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
                raise ValueError("Invalid public F critical values or score_bound")

        s1 = self.engine.sum(x1)
        s2 = self.engine.sum(x2)

        q1 = self.engine.sum(self.engine.mult(x1, x1))
        q2 = self.engine.sum(self.engine.mult(x2, x2))

        s1_squared = self.engine.mult(s1, s1)
        s2_squared = self.engine.mult(s2, s2)

        variance1 = self.engine.sub(q1, self.engine.mult(s1_squared, 1.0 / n1))
        variance2 = self.engine.sub(q2, self.engine.mult(s2_squared, 1.0 / n2))

        variance1 = self.engine.mult(variance1, 1.0 / (n1 - 1))
        variance2 = self.engine.mult(variance2, 1.0 / (n2 - 1))

        var_max = R ** 2 / 4.0

        if var_max <= 0.0:
            raise ValueError("var_max must be positive")

        var2_norm = self.engine.mult(variance2, 1.0 / var_max)
        self._bootstrap(var2_norm)

        inv_sqrt_var2_norm = self.stats.invSqrt(var2_norm)

        inv_var2_norm = self.engine.mult(
            inv_sqrt_var2_norm,
            inv_sqrt_var2_norm,
        )
        inv_var2 = self.engine.mult(inv_var2_norm, 1.0 / var_max)

        f_statistic = self.engine.mult(variance1, inv_var2)

        df1 = n1 - 1
        df2 = n2 - 1

        if lower_critical is None:
            return f_statistic, df1, df2

        inv_sqrt_for_score = Ciphertext(inv_sqrt_var2_norm)
        self._bootstrap(inv_sqrt_for_score)

        inv_var2_for_score = self.engine.mult(
            inv_sqrt_for_score,
            inv_sqrt_for_score,
        )
        inv_var2_for_score = self.engine.mult(
            inv_var2_for_score,
            1.0 / var_max,
        )

        f_for_score = self.engine.mult(variance1, inv_var2_for_score)

        lower_score = self.engine.sub(f_for_score, lower_critical)
        upper_score = self.engine.sub(f_for_score, upper_critical)
        score = self.engine.mult(lower_score, upper_score)

        normalized_score = self.engine.mult(score, 1.0 / score_bound)
        sign_score = self.stats.sign(normalized_score)

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

        if sigma1_sq <= 0.0 or sigma2_sq <= 0.0:
            raise ValueError("sigma1_sq and sigma2_sq must be positive")

        if alpha is None and score_bound is not None:
            raise ValueError("score_bound requires an alpha for decision mode")

        if alpha is not None:
            if not np.isfinite(alpha) or not 0.0 < alpha < 1.0:
                raise ValueError("alpha must be a finite value in (0, 1)")

            if (
                score_bound is None
                or not np.isfinite(score_bound)
                or score_bound <= 0.0
            ):
                raise ValueError(
                    "decision mode requires a finite positive score_bound"
                )

        s1 = self.engine.sum(x1)
        s2 = self.engine.sum(x2)

        mean1 = self.engine.mult(s1, 1.0 / n1)
        mean2 = self.engine.mult(s2, 1.0 / n2)

        numerator = self.engine.sub(mean1, mean2)

        if delta0 != 0.0:
            numerator = self.engine.sub(numerator, delta0)

        v_public = sigma1_sq / n1 + sigma2_sq / n2

        if v_public <= 0.0:
            raise ValueError("V_public must be positive")

        z = self.engine.mult(numerator, 1.0 / np.sqrt(v_public))

        if alpha is None:
            return z

        critical_value = NormalDist().inv_cdf(1.0 - alpha / 2.0)

        z_squared = self.engine.mult(z, z)
        score = self.engine.sub(z_squared, critical_value ** 2)

        normalized_score = self.engine.mult(score, 1.0 / score_bound)
        sign_score = self.stats.sign(normalized_score)

        step = self.engine.add(sign_score, 1.0)
        step = self.engine.mult(step, 0.5)

        return step