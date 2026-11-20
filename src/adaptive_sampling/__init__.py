"""Adaptive sampling / Próbkowanie adaptacyjne. Authors: Andrzej Janaszek, Piotr Wróbel, Maciej Dems, Omer Ceylan, Yasar Gurbuz, Łukasz Kubiszyn, Józef Piotrowski, Rafał Kotyński; code: ChatGPT (OpenAI)."""
from .sampling import SamplingResult, sample_adaptive, interpolate_samples, interval_disagreement
__all__ = ['SamplingResult', 'sample_adaptive', 'interpolate_samples', 'interval_disagreement']
