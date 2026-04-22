"""Device selection and backend availability checks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

Device = Literal["cpu", "cuda", "mps", "auto"]


@dataclass(frozen=True)
class BackendInfo:
    name: str
    device: str
    torch_available: bool
    numba_available: bool


def _try_torch():
    try:
        import torch
        return torch
    except ImportError:
        return None


def _try_numba():
    try:
        import numba
        return numba
    except ImportError:
        return None


def available_devices() -> list[str]:
    devs = ["cpu"]
    torch = _try_torch()
    if torch is not None:
        if torch.cuda.is_available():
            devs.append("cuda")
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            devs.append("mps")
    return devs


def resolve_device(device: Device = "auto") -> str:
    if device != "auto":
        return device
    devs = available_devices()
    for preferred in ("cuda", "mps", "cpu"):
        if preferred in devs:
            return preferred
    return "cpu"


def get_backend(device: Device = "auto") -> BackendInfo:
    resolved = resolve_device(device)
    return BackendInfo(
        name="torch" if resolved != "cpu" or _try_torch() is not None else "numpy",
        device=resolved,
        torch_available=_try_torch() is not None,
        numba_available=_try_numba() is not None,
    )
