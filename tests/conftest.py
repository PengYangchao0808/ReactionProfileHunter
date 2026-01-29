"""
Pytest configuration for ReactionProfileHunter.

This conftest.py ensures rph_core can be imported even when not
installed with pip install -e.
"""
import sys
from pathlib import Path

# Add project root to sys.path for imports (needed when running tests without pip install -e)
project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
