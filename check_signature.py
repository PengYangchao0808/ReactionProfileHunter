import inspect
from rph_core.utils.qc_interface import run_gaussian_optimization

sig = inspect.signature(run_gaussian_optimization)
print(f"def run_gaussian_optimization{sig}")

