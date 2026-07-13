"""V4 canonical PEB scanner."""

from rph_core.steps.step2_retro.peb_engine import PEBScanEngine


class PEBScanner(PEBScanEngine):
    """Canonical S2 engine; only the PEB/path scan entrypoint is supported."""

    def run(self, product_xyz, output_dir, forming_bonds, scan_config=None):
        return super().run(product_xyz, output_dir, forming_bonds, scan_config=scan_config)
