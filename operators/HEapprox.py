import math

import numpy as np
import heaan as hn
from numpy.polynomial import Chebyshev

from engine.HEengine import HEengine
from engine.HEdata import Ciphertext

from coeffs.invSqrt_coeffs import (
    _INV_SQRT_COEFFS,
    _INV_SQRT_MIN,
    _INV_SQRT_MAX,
)
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
        self.engine.bootstrap_extended(ctxt)
        self._bootstrap_count += 1

    def _bootstrap_bounded(
        self,
        ctxt: Ciphertext,
        bound: float,
    ):
        # Bootstrap, performed in place.
        # If bound is greater than 1, perform bootstrap_extended.
        # Else, perform regular bootstrap.
        self.engine.bootstrap_bounded(ctxt, bound)
        self._bootstrap_count += 1

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
        log_degree: int
    ):
        degree = 2**log_degree - 1

        try:
            coeffs = _INV_SQRT_COEFFS[degree]
        except KeyError:
            raise ValueError(
                f"No invSqrt coefficients are available for degree {degree}"
            ) from None

        return np.asarray(coeffs, dtype=np.float64).copy()

    def invSqrt(
        self,
        x_half: Ciphertext,
        x_cheb: Ciphertext,
        log_degree=4,
        iteration=5,
        scale=1.0,
    ):

        # Inputs (prepared by the caller to save levels):
        #
        # x_half = x / 2,
        # x_cheb = 2 * (x - v_min) / (v_max - v_min) - 1,
        # with x in domain = [v_min, v_max].
        y = self.invSqrt_init(x_cheb, log_degree, scale)

        # Refine the initial approximation using Newton iterations.
        for _ in range(iteration):
            y = self.newton_step(x_half, y)

        return y

    def invSqrt_init(
        self,
        x_cheb: Ciphertext,
        log_degree,
        scale
    ):

        # Returns initial approximation y0.
        coeffs = self.inv_sqrt_coeffs(log_degree)

        # Obtain the initial inverse square root approximation
        # from the Chebyshev-mapped ciphertext.
        return self.engine.evaluate_chebyshev(
            x_cheb,
            self.engine.make_cheb_coeffs(coeffs),
            scale
        )

    def newton_step(
        self,
        x_half: Ciphertext,
        y: Ciphertext,
    ):

        if y.level() <= 4:
            self._bootstrap(y)

        tmp_a = self.engine.mult(y, 1.5)

        tmp_b = self.engine.mult(x_half, y)
        y_squared = self.engine.mult(y, y)
        tmp_b = self.engine.mult(tmp_b, y_squared)

        return self.engine.sub(tmp_a, tmp_b)

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
