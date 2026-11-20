"""Adaptive sampling / Próbkowanie adaptacyjne. Author: Rafał Kotyński; code: ChatGPT (OpenAI)."""
from .sampling import SamplingResult, sample_adaptive, interpolate_samples, interval_disagreement
__all__ = ['SamplingResult', 'sample_adaptive', 'interpolate_samples', 'interval_disagreement']
