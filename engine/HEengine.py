import os
import numpy as np
import heaan as hn

from HEdata import Message, Ciphertext


class HEengine:

    def __init__(
        self,
        params=hn.ParameterPreset.FGb,
        device_type="cpu",
        device_id=0,
        log_slots=15,
        setting_root="/root/heaan_setting/",
        separate_keys_by_slots=False,
        warmup_bootstrap=True,
    ):

        self.params = params
        self.device_type = device_type.lower()
        self.device_id = device_id
        self.log_slots = log_slots
        self.num_slots = 2 ** log_slots
        self.separate_keys_by_slots = separate_keys_by_slots
        self.warmup_bootstrap = warmup_bootstrap

        if self.device_type not in ("cpu", "gpu"):
            raise ValueError(
                f"device_type must be 'cpu' or 'gpu': {device_type}"
            )

        params_preset = str(params)[-3:]

        if separate_keys_by_slots:
            self.setting_dir_path = os.path.join(
                setting_root,
                params_preset,
                f"log_slots_{log_slots}",
            )
        else:
            self.setting_dir_path = os.path.join(
                setting_root,
                params_preset,
            )

        self.key_dir_path = os.path.join(
            self.setting_dir_path,
            "keys",
        ) + os.sep

        self.SK_name = "SK"
        self.PK_name = "PK"

        self.context = None
        self.sk = None
        self.pk = None
        self.ect = None
        self.dct = None
        self.evt = None
        self.bts = None
        self.dt = None

        self.__setup__()

    def is_gpu(self):

        return self.device_type == "gpu"

    def __setup__(self):

        if self.is_gpu():
            self.context = hn.make_context(
                self.params,
                {self.device_id},
            )
            self.dt = hn.Device(
                hn.DeviceType.GPU,
                self.device_id,
            )

        else:
            self.context = hn.make_context(self.params)
            self.dt = None

        os.makedirs(self.key_dir_path, exist_ok=True)

        self._load_keys()
        self._init_operators()

        if self.warmup_bootstrap:
            self._warmup_bootstrap()

        return self

    def _load_keys(self):

        sk_path = self.key_dir_path + self.SK_name

        try:
            self.sk = hn.SecretKey(
                self.context,
                sk_path,
            )
            self.pk = hn.KeyPack(
                self.context,
                self.key_dir_path,
            )

        except Exception:
            self.sk = hn.SecretKey(self.context)
            self.sk.save(sk_path)

            keygen = hn.KeyGenerator(
                self.context,
                self.sk,
            )
            keygen.gen_common_keys()
            keygen.gen_rot_keys_for_bootstrap(self.log_slots)
            keygen.save(self.key_dir_path)

            self.pk = keygen.keypack

    def _init_operators(self):

        if self.is_gpu():
            self.sk.to(self.dt)
            self.pk.to(self.dt)

        self.ect = hn.Encryptor(self.context)
        self.dct = hn.Decryptor(self.context)
        self.evt = hn.HomEvaluator(
            self.context,
            self.pk,
        )
        self.bts = hn.Bootstrapper(self.evt)

    def _warmup_bootstrap(self):

        msg = hn.Message(np.zeros(self.num_slots, dtype=np.float64))

        if self.is_gpu():
            msg.to(self.dt)

        ctxt = hn.Ciphertext(self.context)

        self.ect.encrypt(msg, self.pk, ctxt)
        self.bts.bootstrap(ctxt, ctxt)

    def enc(self, msg, level=12):

        if not isinstance(msg, Message):
            raise TypeError("msg must be Message")

        if not isinstance(level, int) or level < 0:
            raise ValueError("level must be a non-negative integer")

        ret = Ciphertext(len(msg), self.context)

        for i in range(len(msg)):
            self.ect.encrypt(msg[i], self.pk, ret[i], level)

        return ret

    def dec(self, ctxt):

        if not isinstance(ctxt, Ciphertext):
            raise TypeError("ctxt must be Ciphertext")

        ret = Message(len(ctxt), self.log_slots)

        for i in range(len(ctxt)):
            self.dct.decrypt(ctxt[i], self.sk, ret[i])
            ret[i].to_host()

        return ret

    def add(self, oprd1, oprd2):

        if not isinstance(oprd1, Ciphertext):
            raise TypeError("oprd1 must be Ciphertext")

        ret = Ciphertext(len(oprd1), self.context)

        if isinstance(oprd2, (Ciphertext, Message)):
            if len(oprd2) == len(oprd1):
                for i in range(len(oprd1)):
                    self.evt.add(oprd1[i], oprd2[i], ret[i])

            elif len(oprd2) == 1:
                for i in range(len(oprd1)):
                    self.evt.add(oprd1[i], oprd2[0], ret[i])

            else:
                raise ValueError(
                    "oprd2 must have length 1 or the same length as oprd1"
                )

        elif isinstance(oprd2, (int, float)):
            for i in range(len(oprd1)):
                self.evt.add(oprd1[i], oprd2, ret[i])

        else:
            raise TypeError(
                "oprd2 must be Ciphertext, Message, int, or float"
            )

        return ret

    def sub(self, oprd1, oprd2):

        if not isinstance(oprd1, Ciphertext):
            raise TypeError("oprd1 must be Ciphertext")

        ret = Ciphertext(len(oprd1), self.context)

        if isinstance(oprd2, (Ciphertext, Message)):
            if len(oprd2) == len(oprd1):
                for i in range(len(oprd1)):
                    self.evt.sub(oprd1[i], oprd2[i], ret[i])

            elif len(oprd2) == 1:
                for i in range(len(oprd1)):
                    self.evt.sub(oprd1[i], oprd2[0], ret[i])

            else:
                raise ValueError(
                    "oprd2 must have length 1 or the same length as oprd1"
                )

        elif isinstance(oprd2, (int, float)):
            for i in range(len(oprd1)):
                self.evt.sub(oprd1[i], oprd2, ret[i])

        else:
            raise TypeError(
                "oprd2 must be Ciphertext, Message, int, or float"
            )

        return ret

    def _mult_integer(self, lhs: Ciphertext, rhs: int):

        if not isinstance(lhs, Ciphertext):
            raise TypeError("lhs must be a Ciphertext")

        if not isinstance(rhs, (int, np.integer)):
            raise TypeError("rhs must be an integer")

        if not hasattr(self.evt, "mult_integer"):
            raise RuntimeError(
                "The current HEaaN evaluator does not support mult_integer"
            )

        ret = Ciphertext(len(lhs), self.context)

        for i in range(len(lhs)):
            self.evt.mult_integer(lhs[i], int(rhs), ret[i])

        return ret

    def mult(self, oprd1, oprd2):

        if not isinstance(oprd1, Ciphertext):
            raise TypeError("oprd1 must be Ciphertext")

        if isinstance(oprd2, (Ciphertext, Message)):
            ret = Ciphertext(len(oprd1), self.context)

            if len(oprd2) == len(oprd1):
                for i in range(len(oprd1)):
                    self.evt.mult(oprd1[i], oprd2[i], ret[i])

            elif len(oprd2) == 1:
                for i in range(len(oprd1)):
                    self.evt.mult(oprd1[i], oprd2[0], ret[i])

            else:
                raise ValueError(
                    "oprd2 must have length 1 or the same length as oprd1"
                )

            return ret

        if isinstance(oprd2, (int, np.integer)):
            return self._mult_integer(oprd1, int(oprd2))

        if isinstance(oprd2, (float, np.floating)):
            ret = Ciphertext(len(oprd1), self.context)

            for i in range(len(oprd1)):
                self.evt.mult(oprd1[i], float(oprd2), ret[i])

            return ret

        raise TypeError(
            "oprd2 must be Ciphertext, Message, int, or float"
        )

    def square(self, oprd):

        if not isinstance(oprd, Ciphertext):
            raise TypeError("oprd must be Ciphertext")

        ret = Ciphertext(len(oprd), self.context)

        for i in range(len(oprd)):
            self.evt.square(oprd[i], ret[i])

        return ret

    def left_rotate(self, ctxt, step):

        if not isinstance(ctxt, Ciphertext):
            raise TypeError("ctxt must be Ciphertext")

        if not isinstance(step, int):
            raise TypeError("step must be int")

        if step < 0:
            raise ValueError("step must be non-negative")

        ret = Ciphertext(len(ctxt), self.context)

        for i in range(len(ctxt)):
            self.evt.left_rotate(ctxt[i], step, ret[i])

        return ret

    def right_rotate(self, ctxt, step):

        if not isinstance(ctxt, Ciphertext):
            raise TypeError("ctxt must be Ciphertext")

        if not isinstance(step, int):
            raise TypeError("step must be int")

        if step < 0:
            raise ValueError("step must be non-negative")

        ret = Ciphertext(len(ctxt), self.context)

        for i in range(len(ctxt)):
            self.evt.right_rotate(ctxt[i], step, ret[i])

        return ret

    def negate(self, ctxt):

        if not isinstance(ctxt, Ciphertext):
            raise TypeError("ctxt must be Ciphertext")

        ret = Ciphertext(len(ctxt), self.context)

        for i in range(len(ctxt)):
            self.evt.negate(ctxt[i], ret[i])

        return ret

    def bootstrap(self, oprd):

        if not isinstance(oprd, Ciphertext):
            raise TypeError("oprd must be Ciphertext")

        for i in range(len(oprd)):
            self.bts.bootstrap(oprd[i], oprd[i])

    def level_down(self, ctxt, target_level):

        if not isinstance(ctxt, Ciphertext):
            raise TypeError("ctxt must be Ciphertext")

        if not isinstance(target_level, int):
            raise TypeError("target_level must be int")

        ret = Ciphertext(ctxt)

        for i in range(len(ret)):
            if ret[i].level > target_level:
                self.evt.level_down(ret[i], target_level, ret[i])

        return ret

    def evaluate_chebyshev(self, ctxt, coeffs, scale=1.0):

        if not isinstance(ctxt, Ciphertext):
            raise TypeError("ctxt must be Ciphertext")

        if not isinstance(coeffs, hn.math.approx.ChebyshevCoefficients):
            raise TypeError("coeffs must be ChebyshevCoefficients instance")

        ret = Ciphertext(len(ctxt), self.context)

        for i in range(len(ctxt)):
            ret.data[i] = hn.math.approx.evaluate_chebyshev_expansion(
                self.evt,
                self.bts,
                ctxt[i],
                coeffs,
                scale,
            )

        return ret

    def sum(self, ctxt: Ciphertext):

        if not isinstance(ctxt, Ciphertext):
            raise TypeError("ctxt must be Ciphertext")

        if len(ctxt) == 0:
            raise ValueError("ctxt must not be empty")

        ret = Ciphertext(1, self.context)
        ret.data[0] = hn.Ciphertext(ctxt[0])

        for i in range(1, len(ctxt)):
            self.evt.add(ret[0], ctxt[i], ret[0])

        for i in range(self.log_slots):
            rotated = self.left_rotate(ret, 2 ** i)
            ret = self.add(ret, rotated)

        return ret