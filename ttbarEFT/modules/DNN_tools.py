import torch

class TensorAccumulator(AccumulatorABC):
    def __init__(self, tensor: torch.Tensor, dtype=torch.float64):
        self._tensor = tensor
        self._dtype = dtype

    def add(self, other: "TensorAccumulator") -> "TensorAccumulator":
        return TensorAccumulator(torch.concat([self._tensor, other._tensor]))

    def __add__(self, other):
        return self.add(other)

    def __iadd__(self, other):
        return self.add(other)

    def get(self) -> torch.Tensor:
        return self._tensor

    def identity(self):
        return TensorAccumulator(torch.Tensor)

    def concat(self, tensor: torch.Tensor):
        return TensorAccumulator(torch.concat([self._tensor, tensor], axis=0).type(self._dtype))
