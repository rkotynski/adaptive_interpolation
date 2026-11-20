"""Adaptive sampling / Próbkowanie adaptacyjne. Author: Rafał Kotyński; code: ChatGPT (OpenAI)."""
from .v7 import V7Result, sample_v7, interpolate_samples, interval_disagreement
__all__ = ['V7Result', 'sample_v7', 'interpolate_samples', 'interval_disagreement']
