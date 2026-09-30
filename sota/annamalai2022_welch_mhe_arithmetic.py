"""Arithmetic trace of Annamalai, Jin, and Aung's 2022 Welch protocol.

Reference: *Communication-Efficient Secure Federated Statistical Tests from
Multiparty Homomorphic Encryption*, Applied Sciences 12(22), 11462 (2022),
doi:10.3390/app122211462.

This implements the paper's encrypted arithmetic: vertically partitioned
``values`` and binary ``labels``, its public ``[0, 1]`` input domain, the paper
upper bounds for the two Welch denominators, and direct (unnormalized)
Goldschmidt HDiv.  HEaaN in this repository has a single secret key, so it
cannot reproduce the paper's MHE collective bootstrap, collective decryption,
or MHE-to-SMPC p-value lookup.  It is therefore an arithmetic trace, not a
cryptographic reimplementation of the full MHE protocol.
"""

from pp_test.engine.HEengine import HEEngine
from hedata.data import HEData
from operators.operator import HEOperator


class Annamalai2022WelchArithmetic:
    """Paper-form Welch ``T^2`` and ``df`` using direct Goldschmidt HDiv."""

    paper = "Annamalai et al. (2022), Algorithm 2 and Section 4.1"
    scope = "single-key HEaaN arithmetic trace; excludes MHE/MPC protocols"

    def __init__(self, engine: HEEngine, verbose: bool = True):
        self._engine = engine
        self._ho = HEOperator(engine)
        self._verbose = verbose

    def bootstrap_count(self) -> int:
        return self._ho.bootstrap_count()

    def reset_bootstrap_count(self) -> None:
        self._ho.reset_bootstrap_count()

    def paper_hdiv(
        self,
        numerator: HEData,
        denominator: HEData,
        denominator_upper: float,
        num_iter: int,
        bootstrap_level: int = 11,
        bootstrap_threshold: int = 4,
    ) -> HEData:
        """Algorithm-2 Goldschmidt division with the paper's ``b_l = 0`` bound.

        For Section 4.1's normalized inputs, the paper gives ``b in [0,b_u]``.
        Substituting ``b_l=0`` in its published linear initializer yields

        ``w_0 = -8*b/b_u^2 + 8/b_u``.

        The valid-input precondition remains ``b > 0`` slotwise: Welch's test
        itself is undefined for zero within-group variation.  No oracle lower
        bound or quotient normalization is introduced here.
        """

        if denominator_upper <= 0.0:
            raise ValueError("denominator_upper must be positive.")
        if num_iter < 1:
            raise ValueError("num_iter must be at least one.")

        ho = self._ho
        w = ho.add_const(
            ho.mult_const(denominator, -8.0 / denominator_upper**2),
            8.0 / denominator_upper,
        )
        a_current = ho.copy_new(numerator)
        b_current = ho.copy_new(denominator)

        for iteration in range(num_iter):
            if self._verbose:
                print(
                    "[Annamalai2022 HDiv] "
                    f"iter={iteration}, a_level={a_current.level()}, "
                    f"b_level={b_current.level()}, w_level={w.level()}"
                )

            # These two HMul operations are independent: depth one per round.
            a_current = ho.mult(a_current, w)
            b_current = ho.mult(b_current, w)
            w = ho.add_const(ho.mult_const(b_current, -1.0), 2.0)

            # Algorithm 2 collectively bootstraps all three live states.  Keep
            # that behavior rather than applying the optimized residual form.
            if w.level() < bootstrap_threshold:
                if self._verbose:
                    print("[Annamalai2022 HDiv] bootstrapping a, b, w")
                a_current = ho.do_bootstrapping(a_current, bootstrap_level)
                b_current = ho.do_bootstrapping(b_current, bootstrap_level)
                w = ho.do_bootstrapping(w, bootstrap_level)

        return a_current

    def welch_t2_df(
        self,
        values: HEData,
        labels: HEData,
        n_label_one: int,
        n_label_zero: int,
        num_iter: int,
        return_debug: bool = False,
    ):
        """Section-4.1 Welch arithmetic from encrypted values and labels.

        ``values`` holds normalized values in ``[0,1]``; ``labels`` holds the
        encrypted binary class vector.  Group sizes are public as required by
        the paper.  The two divisions deliberately use the paper's public
        upper bounds, not data-dependent bounds.
        """

        if values.size() != labels.size():
            raise ValueError("values and labels must have the same length.")
        if n_label_one <= 1 or n_label_zero <= 1:
            raise ValueError("Both public group sizes must exceed one.")
        if n_label_one + n_label_zero != values.size():
            raise ValueError("Public group sizes must sum to values.size().")

        ho = self._ho
        if self._verbose:
            print(f"n(label=1), n(label=0): {n_label_one} {n_label_zero}")

        # Vertical partition: the value owner provides values and the class
        # owner provides labels.  Only encrypted products join the two.
        total_s = ho.sum(values)
        values_squared = ho.mult(values, values)
        total_q = ho.sum(values_squared)
        s_one = ho.sum(ho.mult(values, labels))
        q_one = ho.sum(ho.mult(values_squared, labels))
        s_zero = ho.sub(total_s, s_one)
        q_zero = ho.sub(total_q, q_one)

        mean_one = ho.mult_const(s_one, 1.0 / n_label_one)
        mean_zero = ho.mult_const(s_zero, 1.0 / n_label_zero)
        variance_one = ho.mult_const(
            ho.sub(q_one, ho.mult_const(ho.mult(s_one, s_one), 1.0 / n_label_one)),
            1.0 / (n_label_one - 1),
        )
        variance_zero = ho.mult_const(
            ho.sub(
                q_zero,
                ho.mult_const(ho.mult(s_zero, s_zero), 1.0 / n_label_zero),
            ),
            1.0 / (n_label_zero - 1),
        )

        a_one = ho.mult_const(variance_one, 1.0 / n_label_one)
        a_zero = ho.mult_const(variance_zero, 1.0 / n_label_zero)
        v = ho.add(a_one, a_zero)
        d = ho.add(
            ho.mult_const(ho.mult(a_one, a_one), 1.0 / (n_label_one - 1)),
            ho.mult_const(ho.mult(a_zero, a_zero), 1.0 / (n_label_zero - 1)),
        )
        mean_difference = ho.sub(mean_one, mean_zero)
        t2_numerator = ho.mult(mean_difference, mean_difference)
        df_numerator = ho.mult(v, v)

        # Exact Section-4.1 upper bounds for values in [0, 1].
        v_upper = 1.0 / n_label_one + 1.0 / n_label_zero
        d_upper = (
            1.0 / (n_label_one**2 * (n_label_one - 1))
            + 1.0 / (n_label_zero**2 * (n_label_zero - 1))
        )
        if self._verbose:
            print(f"[Annamalai2022] V upper bound: {v_upper}")
            print(f"[Annamalai2022] D upper bound: {d_upper}")

        t2 = self.paper_hdiv(t2_numerator, v, v_upper, num_iter)
        df = self.paper_hdiv(df_numerator, d, d_upper, num_iter)

        if return_debug:
            return {
                "T2": t2,
                "DF": df,
                "S_one": s_one,
                "S_zero": s_zero,
                "Q_one": q_one,
                "Q_zero": q_zero,
                "V": v,
                "D": d,
                "V_upper": v_upper,
                "D_upper": d_upper,
                "scope": self.scope,
            }
        return t2, df
