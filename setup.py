"""
ReactionProfileHunter v5 - Setup Script
========================================

Author: QCcalc Team
Date: 2026-01-09
"""

from pathlib import Path
from setuptools import setup, find_packages

def read_version() -> str:
    version_file = Path(__file__).parent / "rph_core" / "version.py"
    version_globals = {}
    version_text = version_file.read_text(encoding="utf-8")
    exec(version_text, version_globals)
    return version_globals["__version__"]


# 读取 README
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

setup(
    name="reaction-profile-hunter",
    version=read_version(),
    author="QCcalc Team",
    author_email="your-email@example.com",
    description="过渡态搜索与特征提取的串行工作流系统",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-repo/reaction-profile-hunter",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "rdkit>=2022.09.1",
        "pyyaml>=6.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "black>=22.0.0",
            "mypy>=0.950",
        ],
    },
    entry_points={
        "console_scripts": [
            "rph_run=rph_core.orchestrator:main",
        ],
    },
    include_package_data=True,
    package_data={
        "ReactionProfileHunter": [
            "config/*.yaml",
            "config/templates/*",
        ],
    },
)
