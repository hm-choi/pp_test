import math

import numpy as np
import heaan as hn
from numpy.polynomial import Chebyshev

from engine.HEengine import HEengine
from engine.HEdata import Ciphertext

from coeffs.sign_coeffs import _SIGN_PPTEST_DATA


class HEApprox:

    def __init__(
        self,
        engine: HEengine,
    ):

        self.engine = engine
        self._bootstrap_count = 0

        # Runtime-generated inverse square root coefficients,
        # keyed by (domain, log_degree).
        self._coeff_cache = {}

    def bootstrap_count(
        self,
    ):
        return self._bootstrap_count

    def reset_bootstrap_count(
        self,
    ):

        self._bootstrap_count = 0

    def _bootstrap(
        self,
        ctxt: Ciphertext,
    ):

        # Regular bootstrap, performed in place.
        # Input range [-1, 1], input level >= 3.
        self.engine.bootstrap(ctxt)
        self._bootstrap_count += 1

    def _bootstrap_extended(
        self,
        ctxt: Ciphertext,
    ):

        # Extended bootstrap, performed in place.
        # Input range [-2^20, 2^20], input level >= 4.
        for i in range(len(ctxt)):
            self.engine.bts.bootstrap_extended(ctxt[i], ctxt[i])

        self._bootstrap_count += 1

    def _bootstrap_bounded(
        self,
        ctxt: Ciphertext,
        bound: float,
    ):

        # Choose the bootstrap by the public bound of |ctxt|.
        if bound is not None and bound <= 1.0:
            self._bootstrap(ctxt)
            return

        if ctxt.level() < 4:
            raise ValueError(
                "Extended bootstrap requires level >= 4 "
                f"(got {ctxt.level()})"
            )

        self._bootstrap_extended(ctxt)

    def ensure_level(
        self,
        ctxt: Ciphertext,
        level: int,
        bound=None,
    ):

        # Bootstrap in place if ctxt.level() < level; bound is the public
        # bound of |ctxt| (regular bootstrap if <= 1, extended otherwise).
        if ctxt.level() < level:
            self._bootstrap_bounded(ctxt, bound)

        return ctxt

    def _ensure_cheb_level(
        self,
        ctxt: Ciphertext,
        degree: int,
    ):

        # HEaaN evaluates a degree-(2^k - 1) Chebyshev expansion with
        # k levels and requires input level >= k + 3.
        # The bootstrap is performed in place (input must lie in [-1, 1]).
        log_degree = math.ceil(math.log2(degree + 1))

        if ctxt.level() < log_degree + 3:
            self._bootstrap(ctxt)

    # ------------------------------------------------------------------
    # Inverse square root
    # ------------------------------------------------------------------

    def inv_sqrt_coeffs(
        self,
        log_degree: int,
        domain,
    ):

        # Chebyshev coefficients of 1 / sqrt(x) for x in domain = [v_min, v_max],
        # generated from the public bounds (same construction as
        # coeffs/approx.ipynb):
        # f(t) = 1 / sqrt(x) with x = (v_max - v_min) t / 2 + (v_max + v_min) / 2.
        degree = 2 ** log_degree - 1
        v_min, v_max = float(domain[0]), float(domain[1])

        if not 0.0 < v_min < v_max:
            raise ValueError("domain must satisfy 0 < v_min < v_max")

        key = (v_min, v_max, log_degree)

        if key not in self._coeff_cache:
            def func(t):
                x = (v_max - v_min) * t / 2.0 + (v_max + v_min) / 2.0
                return x ** -0.5

            poly = Chebyshev.interpolate(func, degree, domain=[-1.0, 1.0])
            self._coeff_cache[key] = poly.coef.astype(np.float64)

        return self._coeff_cache[key]

    def invSqrt(
        self,
        x_half: Ciphertext,
        x_cheb: Ciphertext,
        domain,
        log_degree=6,
        iteration=7,
        pre_bts=False,
    ):

        # Inputs (prepared by the caller to save levels):
        #
        # x_half = x / 2,
        # x_cheb = 2x / (v_max - v_min) - (v_max + v_min) / (v_max - v_min),
        #
        # with x in domain = [v_min, v_max].
        # pre_bts: Pre-BTS indicator c of HE-DAP (see invSqrt_init).
        x_half, y = self.invSqrt_init(
            x_half,
            x_cheb,
            domain,
            log_degree,
            pre_bts,
        )

        # Refine the initial approximation using Newton iterations.
        for _ in range(iteration):
            y = self.newton_step(x_half, y)

        return y

    def invSqrt_init(
        self,
        x_half: Ciphertext,
        x_cheb: Ciphertext,
        domain,
        log_degree=6,
        pre_bts=False,
    ):

        # Returns (x_half ready for Newton, initial approximation y0).
        #
        # pre_bts: bootstrap x_cheb (in [-1, 1]) first and derive
        # x_half = ((v_max - v_min) * x_cheb + (v_max + v_min)) / 4 from it,
        # i.e. one bootstrap of the input as in HE-DAP.
        coeffs = self.inv_sqrt_coeffs(log_degree, domain)
        v_min, v_max = float(domain[0]), float(domain[1])

        x_cheb = Ciphertext(x_cheb)

        if pre_bts:
            self._bootstrap(x_cheb)

            x_half = self.engine.mult(x_cheb, (v_max - v_min) / 4.0)
            x_half = self.engine.add(x_half, (v_max + v_min) / 4.0)

        # Obtain the initial inverse square root approximation
        # from the Chebyshev-mapped ciphertext.
        y = self.cheb_invSqrt(x_cheb, coeffs)

        # Public bound of x_half = x / 2 is v_max / 2.
        return self.newton_prepare(x_half, v_max / 2.0), y

    def newton_prepare(
        self,
        x_half: Ciphertext,
        half_bound: float,
    ):

        # The Newton output level is min(level(x_half), level(y)) - 2 and
        # must stay >= 4 for the extended bootstrap of y.
        x_half = Ciphertext(x_half)

        if x_half.level() <= 5:
            self._bootstrap_bounded(x_half, half_bound)

        return x_half

    def newton_step(
        self,
        x_half: Ciphertext,
        y: Ciphertext,
    ):

        # y lies in [1/sqrt(v_max), 1/sqrt(v_min)], outside the [-1, 1]
        # range of the regular bootstrap, so the extended bootstrap is used.
        if y.level() <= 5:
            y = Ciphertext(y)
            self._bootstrap_extended(y)

        tmp_a = self.engine.mult(y, 1.5)

        tmp_b = self.engine.mult(x_half, y)
        y_squared = self.engine.mult(y, y)
        tmp_b = self.engine.mult(tmp_b, y_squared)

        return self.engine.sub(tmp_a, tmp_b)

    def cheb_invSqrt(
        self,
        ctxt: Ciphertext,
        coeffs,
    ):

        x = Ciphertext(ctxt)

        coeffs = np.asarray(coeffs, dtype=np.float64)
        degree = len(coeffs) - 1
        log_degree = math.ceil(math.log2(degree + 1))

        # Ensure that the ciphertext has enough remaining levels
        # before evaluating the Chebyshev polynomial.
        self._ensure_cheb_level(x, degree)

        y = self.engine.evaluate_chebyshev(
            x,
            self.engine._make_cheb_coeffs(coeffs),
        )

        # Refresh y0 with the regular bootstrap when the first Newton step
        # would need a bootstrap. y0 may exceed [-1, 1] (up to 1/sqrt(v_min));
        # the regular bootstrap then adds a relative error of about
        # 1e-7 * |y0|^2 (1e-4 at |y0| = 31.6), which the following Newton
        # iterations refine, so the extended bootstrap is not needed here.
        if y.level() <= 5:
            self._bootstrap(y)

        return y

    # ------------------------------------------------------------------
    # Sign
    # ------------------------------------------------------------------

    def sign(
        self,
        ctxt: Ciphertext,
    ):

        # Composite sign approximation; input must lie in [-1, 1].
        ret = Ciphertext(ctxt)

        for coeffs in _SIGN_PPTEST_DATA:
            # Ensure that enough ciphertext levels remain before every
            # polynomial evaluation (the output of each stage stays in
            # [-1, 1], so the regular bootstrap can be used).
            degree = len(coeffs) - 1
            self._ensure_cheb_level(ret, degree)

            ret = self.engine.evaluate_chebyshev(
                ret,
                self.engine._make_cheb_coeffs(coeffs),
            )

        return ret
