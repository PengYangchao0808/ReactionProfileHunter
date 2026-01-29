"""
测试沙盒执行机制 - 验证 XTB/CREST 在有毒路径下的运行
===========================================

该测试验证在包含空格和特殊字符的路径下，QC 计算能够正常进行。

Author: QC Descriptors Team
Date: 2026-01-13
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from rph_core.utils.qc_interface import is_path_toxic, run_in_sandbox
from rph_core.utils.qc_interface import XTBInterface, CRESTInterface
import sys

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class TestToxicPathDetection:
    """测试有毒路径检测"""

    def test_path_with_space(self):
        """测试路径包含空格"""
        path = Path("/home/user/my calculations/test")
        assert is_path_toxic(path) is True

    def test_path_with_brackets(self):
        """测试路径包含方括号"""
        path = Path("/home/user/[test]/results")
        assert is_path_toxic(path) is True

    def test_clean_path(self):
        """测试干净路径"""
        path = Path("/home/user/calc/test")
        assert is_path_toxic(path) is False

    def test_toxic_chars_detection(self):
        """测试各种有毒字符"""
        toxic_chars = [' ', '[', ']', '(', ')', '{', '}']
        for char in toxic_chars:
            path = Path(f"/home/user/test{char}dir")
            assert is_path_toxic(path) is True


class TestSandboxExecution:
    """测试沙盒执行机制"""

    @pytest.fixture
    def sample_xyz(self):
        """创建简单的 XYZ 文件用于测试"""
        return """3
Water
O  0.0000  0.0000  0.0000
H  0.7586  0.0000  0.5044
H -0.7586  0.0000  0.5044
"""

    @pytest.fixture
    def sample_xyz_file(self, sample_xyz, tmp_path):
        """创建测试用的 XYZ 文件"""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text(sample_xyz)
        return xyz_file

    def test_sandbox_basic_execution(self, sample_xyz_file, tmp_path):
        """测试基本的沙盒执行"""
        temp_dir = tmp_path / "sandbox_test"
        temp_dir.mkdir()

        result = run_in_sandbox(
            cmd=["echo", "test"],
            input_files={"input.xyz": sample_xyz_file},
            output_files=["output.xyz"],
            output_dir=temp_dir,
            timeout=10
        )

        assert result.returncode == 0
        assert (temp_dir / "input.xyz").exists()
        print("✓ 基本沙盒执行测试通过")

    def test_sandbox_with_toxic_path(self, sample_xyz_file):
        """测试在有毒路径下使用沙盒"""
        with tempfile.TemporaryDirectory(prefix='test_toxic_') as base_dir:
            toxic_base = Path(base_dir) / "[5+2] Test Space"
            toxic_base.mkdir()
            output_dir = toxic_base / "output"
            output_dir.mkdir()

            result = run_in_sandbox(
                cmd=["echo", "test"],
                input_files={"input.xyz": sample_xyz_file},
                output_files=["output.xyz"],
                output_dir=output_dir,
                timeout=10
            )

            assert result.returncode == 0
            print("✓ 有毒路径沙盒执行测试通过")


class TestXTBWithSandbox:
    """测试 XTB 在沙盒中的运行"""

    @pytest.fixture
    def xtb_interface(self):
        """初始化 XTB 接口"""
        return XTBInterface(
            gfn_level=2,
            solvent="water",
            nproc=2
        )

    @pytest.fixture
    def sample_xyz(self):
        """创建简单的 XYZ 文件用于测试"""
        return """3
Water
O  0.0000  0.0000  0.0000
H  0.7586  0.0000  0.5044
H -0.7586  0.0000  0.5044
"""

    @pytest.fixture
    def sample_xyz_file(self, sample_xyz, tmp_path):
        """创建测试用的 XYZ 文件"""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text(sample_xyz)
        return xyz_file

    def test_xtb_optimize_with_toxic_path(self, xtb_interface, sample_xyz_file):
        """测试 XTB 优化在有毒路径下运行（使用沙盒）"""
        pytest.skip("需要 XTB 可执行文件")

        with tempfile.TemporaryDirectory(prefix='test_xtb_toxic_') as base_dir:
            toxic_base = Path(base_dir) / "[5+2] XTB Test"
            toxic_base.mkdir()
            output_dir = toxic_base / "output"

            result = xtb_interface.optimize(
                xyz_file=sample_xyz_file,
                output_dir=output_dir
            )

            assert result.converged
            assert (output_dir / "xtbopt.xyz").exists()
            print(f"✓ XTB 在有毒路径下优化成功: {output_dir}")


class TestCRESTWithSandbox:
    """测试 CREST 在沙盒中的运行"""

    @pytest.fixture
    def crest_interface(self):
        """初始化 CREST 接口"""
        return CRESTInterface(
            gfn_level=2,
            solvent="water",
            nproc=2
        )

    @pytest.fixture
    def sample_xyz(self):
        """创建简单的 XYZ 文件用于测试"""
        return """3
Water
O  0.0000  0.0000  0.0000
H  0.7586  0.0000  0.5044
H -0.7586  0.0000  0.5044
"""

    @pytest.fixture
    def sample_xyz_file(self, sample_xyz, tmp_path):
        """创建测试用的 XYZ 文件"""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text(sample_xyz)
        return xyz_file

    def test_crest_with_toxic_path(self, crest_interface, sample_xyz_file):
        """测试 CREST 在有毒路径下运行（使用沙盒）"""
        pytest.skip("需要 CREST 可执行文件")

        with tempfile.TemporaryDirectory(prefix='test_crest_toxic_') as base_dir:
            toxic_base = Path(base_dir) / "[5+2] CREST Test Space"
            toxic_base.mkdir()
            output_dir = toxic_base / "output"

            result = crest_interface.run_conformer_search(
                xyz_file=sample_xyz_file,
                output_dir=output_dir
            )

            assert result.exists()
            print(f"✓ CREST 在有毒路径下运行成功: {output_dir}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
