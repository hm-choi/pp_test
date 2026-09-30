import math

import numpy as np
import heaan as hn

from engine.HEengine import HEengine
from engine.HEdata import Ciphertext
from engine.chebyshev_eval import next_power_of_two_sqrt, make_hn_cheb

from coeffs.invSqrt_coeffs import _INV_SQRT_COEFFS
from coeffs.sqrt_coeffs import _SQRT_COEFFS
from coeffs.sign_coeffs import _SIGN_PPTEST_DATA


class HEStats:

    def __init__(self, engine: HEengine):

        self.engine = engine
        self._bootstrap_count = 0

    def bootstrap_count(self):

        return self._bootstrap_count

    def reset_bootstrap_count(self):

        self._bootstrap_count = 0

    def _bootstrap(self, ctxt: Ciphertext):

        self.engine.bootstrap(ctxt)
        self._bootstrap_count += 1

    def _make_cheb_coeffs(
        self,
        coeffs,
    ) -> hn.math.approx.ChebyshevCoefficients:
        coeffs = np.asarray(coeffs, dtype=np.float64)

        degree = len(coeffs) - 1
        baby_step = next_power_of_two_sqrt(degree)

        return make_hn_cheb(coeffs, baby_step)

    def _load_coeffs(self, prefix: str, log_degree: int):

        degree = 2 ** log_degree - 1

        if prefix == "inv":
            coeffs_data = _INV_SQRT_COEFFS
        elif prefix == "sqrt":
            coeffs_data = _SQRT_COEFFS
        else:
            raise ValueError(f"Unsupported coefficient prefix: {prefix}")

        if degree not in coeffs_data:
            raise ValueError(f"Unsupported degree: {degree}")

        return self._make_cheb_coeffs(coeffs_data[degree])

    def invSqrt(
        self,
        ctxt: Ciphertext,
        log_degree=6,
        input_scale=1.0,
        output_scale=1.0,
        iteration=7,
    ):

        coeffs = self._load_coeffs("inv", log_degree)

        x = Ciphertext(ctxt)

        if input_scale != 1.0:
            x = self.engine.mult(x, input_scale)

        y = self.cheb_invSqrt(x, coeffs, output_scale)

        self._bootstrap(y)

        return self.he_newton(x, y, iteration)

    def he_newton(
        self,
        x: Ciphertext,
        y: Ciphertext,
        iteration=10,
    ):

        # Newton iteration for inverse square root:
        # y <- 1.5y - 0.5xy^3
        x = self.engine.mult(x, 0.5)

        if x.level() <= 4:
            self._bootstrap(x)

        for _ in range(iteration):
            if y.level() <= 4:
                self._bootstrap(y)

            tmp_a = self.engine.mult(y, 1.5)

            tmp_b = self.engine.mult(x, y)
            y_sqr = self.engine.mult(y, y)
            tmp_b = self.engine.mult(tmp_b, y_sqr)

            y = self.engine.sub(tmp_a, tmp_b)

        return y

    def cheb_invSqrt(
        self,
        ctxt: Ciphertext,
        coeffs: hn.math.approx.ChebyshevCoefficients,
        output_scale=1.0,
    ):

        x = Ciphertext(ctxt)

        # Shift the input to the Chebyshev approximation domain.
        x = self.engine.sub(x, 1.0)

        ret = self.engine.evaluate_chebyshev(x, coeffs)

        if output_scale != 1.0:
            ret = self.engine.mult(ret, output_scale)

        return ret

    def sqrt(
        self,
        ctxt: Ciphertext,
        log_degree=8,
        input_scale=1.0,
        output_scale=1.0,
    ):

        coeffs = self._load_coeffs("sqrt", log_degree)

        return self.cheb_sqrt(ctxt, coeffs, input_scale, output_scale)

    def cheb_sqrt(
        self,
        ctxt: Ciphertext,
        coeffs: hn.math.approx.ChebyshevCoefficients,
        input_scale=1.0,
        output_scale=1.0,
    ):

        x = Ciphertext(ctxt)

        if input_scale != 1.0:
            x = self.engine.mult(x, input_scale)

        x = self.engine.sub(x, 1.0)

        ret = self.engine.evaluate_chebyshev(x, coeffs)

        if output_scale != 1.0:
            ret = self.engine.mult(ret, output_scale)

        return ret

    def sign(self, ctxt: Ciphertext):

        coeffs_list = [
            self._make_cheb_coeffs(coeffs)
            for coeffs in _SIGN_PPTEST_DATA
        ]

        ret = Ciphertext(ctxt)

        for i, coeffs in enumerate(coeffs_list):
            ret = self.engine.evaluate_chebyshev(ret, coeffs)

            if i != len(coeffs_list) - 1:
                next_degree = len(_SIGN_PPTEST_DATA[i + 1]) - 1
                next_depth = math.ceil(math.log2(next_degree))

                if ret.level() - next_depth < 3:
                    self._bootstrap(ret)

        return ret