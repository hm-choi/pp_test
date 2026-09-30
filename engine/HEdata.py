import numpy as np
import heaan as hn


class Message:

    def __init__(self, input, log_slots=15, value=0):

        num_slots = 2 ** log_slots

        if isinstance(input, np.ndarray):
            self.data = []
            for i in range(0, len(input), num_slots):
                chunk = input[i:i + num_slots]
                if len(chunk) < num_slots:
                    padded = np.zeros(num_slots, dtype=input.dtype)
                    padded[:len(chunk)] = chunk
                    chunk = padded
                self.data.append(hn.Message(chunk))

        elif isinstance(input, int):
            if not isinstance(value, (int, float)):
                raise TypeError("value must be int or float")
            
            self.data = [hn.Message(np.full(num_slots, value)) for _ in range(input)]

        else:
            raise TypeError("input must be np.ndarray or int")

    def __getitem__(self, idx): return self.data[idx]

    def __len__(self): return len(self.data)

    def to(self, device):

        if not isinstance(device, hn.Device):
            raise TypeError("device must be hn.Device")

        for m in self.data: m.to(device)

        return self


class Ciphertext:

    def __init__(self, input, context=None):
        
        if isinstance(input, int):
            if not isinstance(context, hn.Context):
                raise TypeError("context must be hn.Context")
            
            self.data = [hn.Ciphertext(context) for _ in range(input)]

        elif isinstance(input, Ciphertext):
            self.data = [hn.Ciphertext(c) for c in input.data]

        else:
            raise TypeError("input must be int or Ciphertext")

    def __getitem__(self, idx): return self.data[idx]

    def __len__(self): return len(self.data)

    def to(self, device):

        if not isinstance(device, hn.Device):
            raise TypeError("device must be hn.Device")

        for c in self.data: c.to(device)

        return self

    def level(self):
        
        if not self.data:
            raise ValueError("Ciphertext is empty")

        level = self.data[0].level
        for c in self.data[1:]:
            if c.level != level:
                raise ValueError("all ciphertexts must have the same level")

        return level