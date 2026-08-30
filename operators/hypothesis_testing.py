from engine.engine import HEEngine
import heaan as hn
import numpy as np
from hedata.data import HEData
from operators.inv_sqrt import HEStatistics
from operators.operator import HEOperator
import math, json


class HEHypothesisTesting:
    def __init__(self, engine: HEEngine):
        self._engine = engine 
        self.__ho = HEOperator(engine)
        self.__hs = HEStatistics(engine)

    def HE_Welch_T_Test(self, x1: HEData, x2: HEData, R, return_debug=False):
        """
        HE-friendly Welch's t-test.

        x1, x2:
            각 그룹의 연속형 데이터 ciphertext.
            예: smoker=yes charges/1000, smoker=no charges/1000

        R:
            데이터 range upper bound.
            예: charges/1000 이 [0, 63.8] 근처이면 R=63.8 사용 가능.

        return_debug:
            True이면 T뿐 아니라 mean, var, V 등 중간값도 함께 반환.
        """

        ho = self.__ho
        hs = self.__hs

        n1 = x1.size()
        n2 = x2.size()

        print("n1, n2:", n1, n2)

        if n1 <= 1 or n2 <= 1:
            raise ValueError("Welch t-test requires n1 > 1 and n2 > 1.")

        # ============================================================
        # 1. Sufficient statistics
        #    S = Σx
        #    Q = Σx²
        # ============================================================

        S1 = ho.sum(x1)
        S2 = ho.sum(x2)

        x1_sq = ho.mult(x1, x1)
        x2_sq = ho.mult(x2, x2)

        Q1 = ho.sum(x1_sq)
        Q2 = ho.sum(x2_sq)

        # ============================================================
        # 2. Mean
        #    mean = S / n
        # ============================================================

        mean_x1 = ho.mult_const(S1, 1.0 / n1)
        mean_x2 = ho.mult_const(S2, 1.0 / n2)

        # ============================================================
        # 3. Sample variance
        #
        #    s² = (Q - S²/n) / (n-1)
        #
        #    이 방식은 padding slot에 mean을 빼지 않으므로 안전함.
        # ============================================================

        S1_sq = ho.mult(S1, S1)
        S2_sq = ho.mult(S2, S2)

        S1_sq_over_n = ho.mult_const(S1_sq, 1.0 / n1)
        S2_sq_over_n = ho.mult_const(S2_sq, 1.0 / n2)

        var1_num = ho.sub(Q1, S1_sq_over_n)
        var2_num = ho.sub(Q2, S2_sq_over_n)

        var1 = ho.mult_const(var1_num, 1.0 / (n1 - 1))
        var2 = ho.mult_const(var2_num, 1.0 / (n2 - 1))

        # ============================================================
        # 4. Welch denominator inside sqrt
        #
        #    V = s1²/n1 + s2²/n2
        # ============================================================

        V1 = ho.mult_const(var1, 1.0 / n1)
        V2 = ho.mult_const(var2, 1.0 / n2)

        V = ho.add(V1, V2)

        # ============================================================
        # 5. Public upper bound for V
        #
        #    If x ∈ [0, R], then sample variance ≤ R²/4 approximately.
        #
        #    V_max = R²/4 * (1/(n1-1) + 1/(n2-1))
        #
        #    More conservative than empirical bound.
        # ============================================================

        V_max = ((R ** 2) / 4.0) * (
            1.0 / (n1 - 1) + 1.0 / (n2 - 1)
        )

        print("V_max:", V_max)

        if V_max <= 0:
            raise ValueError("V_max must be positive.")

        # ============================================================
        # 6. Normalize V into (0, 1]
        #
        #    V_norm = V / V_max
        #
        #    InvSqrt(V) = InvSqrt(V_norm) / sqrt(V_max)
        # ============================================================

        V_norm = ho.mult_const(V, 1.0 / V_max)

        # 필요하면 InvSqrt 전에 bootstrapping
        V_norm = ho.do_bootstrapping(V_norm, 11)

        inv_sqrt_V_norm = hs.he_inv_sqrt(V_norm)

        # ============================================================
        # 7. Welch t-statistic
        #
        #    t = (mean1 - mean2) / sqrt(V)
        #      = (mean1 - mean2) * InvSqrt(V_norm) / sqrt(V_max)
        # ============================================================

        numerator = ho.sub(mean_x1, mean_x2)
        numerator = ho.mult_const(numerator, 1.0 / np.sqrt(V_max))

        T = ho.mult(numerator, inv_sqrt_V_norm)

        if return_debug:
            return {
                "T": T,
                "mean_x1": mean_x1,
                "mean_x2": mean_x2,
                "var1": var1,
                "var2": var2,
                "V": V,
                "V_norm": V_norm,
                "V_max": V_max,
                "S1": S1,
                "S2": S2,
                "Q1": Q1,
                "Q2": Q2,
            }

        return T

    def HE_Welch_T_Test_With_DF(self, x1: HEData, x2: HEData, R, return_debug=False):
        """
        HE-friendly Welch's t-test with Welch-Satterthwaite degrees of freedom.

        Returns:
            T, DF
        """

        ho = self.__ho
        hs = self.__hs

        n1 = x1.size()
        n2 = x2.size()

        print("n1, n2:", n1, n2)

        if n1 <= 1 or n2 <= 1:
            raise ValueError("Welch t-test requires n1 > 1 and n2 > 1.")

        # ============================================================
        # 1. Sufficient statistics
        #    S = Σx
        #    Q = Σx²
        # ============================================================

        S1 = ho.sum(x1)
        S2 = ho.sum(x2)

        x1_sq = ho.mult(x1, x1)
        x2_sq = ho.mult(x2, x2)

        Q1 = ho.sum(x1_sq)
        Q2 = ho.sum(x2_sq)

        # ============================================================
        # 2. Mean
        # ============================================================

        mean_x1 = ho.mult_const(S1, 1.0 / n1)
        mean_x2 = ho.mult_const(S2, 1.0 / n2)

        # ============================================================
        # 3. Sample variance
        #
        #    s² = (Q - S²/n) / (n-1)
        # ============================================================

        S1_sq = ho.mult(S1, S1)
        S2_sq = ho.mult(S2, S2)

        S1_sq_over_n = ho.mult_const(S1_sq, 1.0 / n1)
        S2_sq_over_n = ho.mult_const(S2_sq, 1.0 / n2)

        var1_num = ho.sub(Q1, S1_sq_over_n)
        var2_num = ho.sub(Q2, S2_sq_over_n)

        var1 = ho.mult_const(var1_num, 1.0 / (n1 - 1))
        var2 = ho.mult_const(var2_num, 1.0 / (n2 - 1))

        # ============================================================
        # 4. Welch denominator
        #
        #    A = s1²/n1
        #    B = s2²/n2
        #    V = A + B
        # ============================================================

        A = ho.mult_const(var1, 1.0 / n1)
        B = ho.mult_const(var2, 1.0 / n2)

        V = ho.add(A, B)

        # ============================================================
        # 5. Public upper bound for V
        #
        #    V_max = R²/4 * (1/(n1-1) + 1/(n2-1))
        # ============================================================

        V_max = ((R ** 2) / 4.0) * (
            1.0 / (n1 - 1) + 1.0 / (n2 - 1)
        )

        print("V_max:", V_max)

        if V_max <= 0:
            raise ValueError("V_max must be positive.")

        # ============================================================
        # 6. Normalize V and compute T
        #
        #    V_norm = V / V_max
        #    InvSqrt(V) = InvSqrt(V_norm) / sqrt(V_max)
        # ============================================================

        V_norm = ho.mult_const(V, 1.0 / V_max)

        V_norm = ho.do_bootstrapping(V_norm, 11)

        inv_sqrt_V_norm = hs.he_inv_sqrt(V_norm)

        numerator = ho.sub(mean_x1, mean_x2)
        numerator = ho.mult_const(numerator, 1.0 / np.sqrt(V_max))

        T = ho.mult(numerator, inv_sqrt_V_norm)

        # ============================================================
        # 7. Welch-Satterthwaite degrees of freedom
        #
        #    df = V² / D
        #
        #    D = A²/(n1-1) + B²/(n2-1)
        # ============================================================

        A_sq = ho.mult(A, A)
        B_sq = ho.mult(B, B)

        D1 = ho.mult_const(A_sq, 1.0 / (n1 - 1))
        D2 = ho.mult_const(B_sq, 1.0 / (n2 - 1))

        D = ho.add(D1, D2)

        # ============================================================
        # 8. Public upper bound for D
        #
        #    A_max = R²/4 * 1/(n1-1)
        #    B_max = R²/4 * 1/(n2-1)
        #
        #    D_max = A_max²/(n1-1) + B_max²/(n2-1)
        #          = (R^4 / 16) * (1/(n1-1)^3 + 1/(n2-1)^3)
        # ============================================================

        A_max = ((R ** 2) / 4.0) * (1.0 / (n1 - 1))
        B_max = ((R ** 2) / 4.0) * (1.0 / (n2 - 1))

        D_max = (A_max ** 2) / (n1 - 1) + (B_max ** 2) / (n2 - 1)

        print("D_max:", D_max)

        if D_max <= 0:
            raise ValueError("D_max must be positive.")

        # ============================================================
        # 9. Normalize D and compute 1/D
        #
        #    D_norm = D / D_max
        #    1/D = InvSqrt(D_norm)^2 / D_max
        # ============================================================

        D_norm = ho.mult_const(D, 1.0 / D_max)

        D_norm = ho.do_bootstrapping(D_norm, 11)

        inv_sqrt_D_norm = hs.he_inv_sqrt(D_norm)

        inv_D_norm = ho.mult(inv_sqrt_D_norm, inv_sqrt_D_norm)

        # ============================================================
        # 10. Compute df
        #
        #    df = V² / D
        #
        #    V = V_max * V_norm
        #    D = D_max * D_norm
        #
        #    df = V_norm² * Inv(D_norm) * V_max² / D_max
        # ============================================================

        V_norm_sq = ho.mult(V_norm, V_norm)

        DF = ho.mult(V_norm_sq, inv_D_norm)
        DF = ho.mult_const(DF, (V_max ** 2) / D_max)

        if return_debug:
            return {
                "T": T,
                "DF": DF,
                "mean_x1": mean_x1,
                "mean_x2": mean_x2,
                "var1": var1,
                "var2": var2,
                "A": A,
                "B": B,
                "V": V,
                "V_norm": V_norm,
                "D": D,
                "D_norm": D_norm,
                "V_max": V_max,
                "D_max": D_max,
            }

        return T, DF

    def HE_F_Test(self, x1: HEData, x2: HEData, R, return_debug=False):
        """
        HE-friendly two-sample variance-ratio F-test.

        Test:
            H0: sigma1^2 = sigma2^2
            H1: sigma1^2 != sigma2^2

        Statistic:
            F = s1^2 / s2^2

        Degrees of freedom:
            df1 = n1 - 1
            df2 = n2 - 1

        Notes:
            - This computes fixed-orientation F = var1 / var2.
            - Do not compute max(var1,var2)/min(var1,var2) under HE.
            - Two-sided p-value can be computed after decryption.
        """

        ho = self.__ho
        hs = self.__hs

        n1 = x1.size()
        n2 = x2.size()

        print("n1, n2:", n1, n2)

        if n1 <= 1 or n2 <= 1:
            raise ValueError("F-test requires n1 > 1 and n2 > 1.")

        # ============================================================
        # 1. Sufficient statistics
        #    S = Σx
        #    Q = Σx²
        # ============================================================

        S1 = ho.sum(x1)
        S2 = ho.sum(x2)

        x1_sq = ho.mult(x1, x1)
        x2_sq = ho.mult(x2, x2)

        Q1 = ho.sum(x1_sq)
        Q2 = ho.sum(x2_sq)

        # ============================================================
        # 2. Mean, optional but useful for debugging
        # ============================================================

        mean_x1 = ho.mult_const(S1, 1.0 / n1)
        mean_x2 = ho.mult_const(S2, 1.0 / n2)

        # ============================================================
        # 3. Sample variance
        #
        #    s² = (Q - S²/n) / (n - 1)
        # ============================================================

        S1_sq = ho.mult(S1, S1)
        S2_sq = ho.mult(S2, S2)

        S1_sq_over_n = ho.mult_const(S1_sq, 1.0 / n1)
        S2_sq_over_n = ho.mult_const(S2_sq, 1.0 / n2)

        var1_num = ho.sub(Q1, S1_sq_over_n)
        var2_num = ho.sub(Q2, S2_sq_over_n)

        var1 = ho.mult_const(var1_num, 1.0 / (n1 - 1))
        var2 = ho.mult_const(var2_num, 1.0 / (n2 - 1))

        # ============================================================
        # 4. Compute F = var1 / var2
        #
        #    Use InvSqrt(var2)^2 instead of direct division.
        # ============================================================

        # Public upper bound for sample variance.
        # If x in [0, R], then s² <= approximately R²/4.
        var_max = (R ** 2) / 4.0

        print("var_max:", var_max)

        if var_max <= 0:
            raise ValueError("var_max must be positive.")

        # Normalize var2 into (0, 1]
        var2_norm = ho.mult_const(var2, 1.0 / var_max)

        # Optional bootstrapping before inverse square-root
        var2_norm = ho.do_bootstrapping(var2_norm, 11)

        inv_sqrt_var2_norm = hs.he_inv_sqrt(var2_norm)

        # 1/var2 = InvSqrt(var2_norm)^2 / var_max
        inv_var2_norm = ho.mult(inv_sqrt_var2_norm, inv_sqrt_var2_norm)
        inv_var2 = ho.mult_const(inv_var2_norm, 1.0 / var_max)

        F = ho.mult(var1, inv_var2)

        # df는 public scalar이므로 굳이 ciphertext로 만들 필요 없음
        df1 = n1 - 1
        df2 = n2 - 1

        if return_debug:
            return {
                "F": F,
                "df1": df1,
                "df2": df2,
                "mean_x1": mean_x1,
                "mean_x2": mean_x2,
                "var1": var1,
                "var2": var2,
                "var2_norm": var2_norm,
                "var_max": var_max,
                "S1": S1,
                "S2": S2,
                "Q1": Q1,
                "Q2": Q2,
            }

        return F, df1, df2

    def HE_Z_Test(self, x1: HEData, x2: HEData, sigma1_sq: float, sigma2_sq: float, delta0: float = 0.0, return_debug: bool = False):
        """
        HE-friendly two-sample Z-test with known/public variances.

        Hypotheses:
            H0: mu1 - mu2 = delta0
            H1: mu1 - mu2 != delta0

        Statistic:
            Z = (mean1 - mean2 - delta0) / sqrt(sigma1_sq/n1 + sigma2_sq/n2)

        Notes:
            - sigma1_sq and sigma2_sq are public/external variances.
            - This version does not require HE InvSqrt.
            - p-value is computed after decryption using standard normal distribution.
        """

        ho = self.__ho

        n1 = x1.size()
        n2 = x2.size()

        print("n1, n2:", n1, n2)

        if n1 <= 0 or n2 <= 0:
            raise ValueError("Z-test requires n1 > 0 and n2 > 0.")

        if sigma1_sq <= 0 or sigma2_sq <= 0:
            raise ValueError("sigma1_sq and sigma2_sq must be positive.")

        # ============================================================
        # 1. Sufficient statistics
        #    S = Σx
        # ============================================================

        S1 = ho.sum(x1)
        S2 = ho.sum(x2)

        # ============================================================
        # 2. Mean
        #    mean = S / n
        # ============================================================

        mean_x1 = ho.mult_const(S1, 1.0 / n1)
        mean_x2 = ho.mult_const(S2, 1.0 / n2)

        # ============================================================
        # 3. Numerator
        #    mean1 - mean2 - delta0
        # ============================================================

        numerator = ho.sub(mean_x1, mean_x2)

        if delta0 != 0.0:
            numerator = ho.add_const(numerator, -delta0)

        # ============================================================
        # 4. Public denominator
        #
        #    V = sigma1_sq/n1 + sigma2_sq/n2
        # ============================================================

        V_public = sigma1_sq / n1 + sigma2_sq / n2

        if V_public <= 0:
            raise ValueError("V_public must be positive.")

        inv_sqrt_V = 1.0 / np.sqrt(V_public)

        Z = ho.mult_const(numerator, inv_sqrt_V)

        if return_debug:
            return {
                "Z": Z,
                "mean_x1": mean_x1,
                "mean_x2": mean_x2,
                "V_public": V_public,
                "sigma1_sq": sigma1_sq,
                "sigma2_sq": sigma2_sq,
                "delta0": delta0,
            }

        return Z