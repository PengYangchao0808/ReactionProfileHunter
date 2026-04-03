"""
测试 ORCA Interface 模块
==========================

Session #1: 测试 _generate_input() 方法
Session #2: 测试 _parse_output() 方法
Session #3: 测试 _find_orca_binary() + _run_orca() 方法
Session #4: 测试 single_point() 端到端流程

Author: QC Descriptors Team
Date: 2026-01-10
"""

import pytest
from pathlib import Path
from typing import cast
from unittest.mock import Mock, patch, MagicMock, call
from rph_core.utils.orca_interface import ORCAInterface, QCResult


@pytest.fixture
def sample_output(tmp_path):
    """模拟 ORCA 输出文件"""
    out_file = tmp_path / "orca.out"
    out_file.write_text("""
ORCA CALCULATION DONE
FINAL SINGLE POINT ENERGY      -123.45678901
ORCA TERMINATED NORMALLY
""")
    return out_file


@pytest.fixture
def sample_xyz(tmp_path):
    """模拟 XYZ 文件"""
    xyz_file = tmp_path / "test.xyz"
    xyz_file.write_text("""3
comment
C    0.0  0.0  0.0
H    1.0  0.0  0.0
H   -1.0  0.0  0.0
""")
    return xyz_file


def test_generate_input_m062x(sample_xyz, tmp_path):
    """测试 M062X 输入生成"""
    orca = ORCAInterface(method="M062X", basis="def2-TZVPP")
    inp_file = orca._generate_input(sample_xyz, tmp_path)

    content = inp_file.read_text()
    assert "! M062X def2-TZVPP def2/J RIJCOSX" in content
    assert "%maxcore 4000" in content
    assert "%pal nprocs 16" in content


def test_generate_input_with_solvent(sample_xyz, tmp_path):
    """测试溶剂模型输入生成"""
    orca = ORCAInterface(solvent="acetone")
    inp_file = orca._generate_input(sample_xyz, tmp_path)

    content = inp_file.read_text()
    assert "%cpcm" in content
    assert "SMDsolvent" in content
    assert "acetone" in content.lower()
    assert "smd true" in content


def test_generate_input_multiple_solvents(sample_xyz, tmp_path):
    """测试多种溶剂支持"""
    solvents = ["water", "toluene", "dichloromethane", "ethanol"]

    for solvent in solvents:
        orca = ORCAInterface(solvent=solvent)
        inp_file = orca._generate_input(sample_xyz, tmp_path)

        content = inp_file.read_text()
        assert "%cpcm" in content
        assert "SMDsolvent" in content
        assert solvent.lower() in content.lower()
        assert "smd true" in content


def test_generate_input_no_solvent(sample_xyz, tmp_path):
    """测试无溶剂情况"""
    orca = ORCAInterface(solvent=cast(str, cast(object, None)))
    inp_file = orca._generate_input(sample_xyz, tmp_path)

    content = inp_file.read_text()
    # 不应该包含 %cpcm 块
    assert "%cpcm" not in content
    assert "SMDsolvent" not in content


def test_generate_input_solvent_none_string(sample_xyz, tmp_path):
    """测试 solvent='None' 字符串"""
    orca = ORCAInterface(solvent="None")
    inp_file = orca._generate_input(sample_xyz, tmp_path)

    content = inp_file.read_text()
    # 不应该包含 %cpcm 块（字符串"None"被转换为None处理）
    assert "%cpcm" not in content


def test_generate_input_double_hybrid(sample_xyz, tmp_path):
    """测试双杂化泛函自动添加 /C 基组"""
    orca = ORCAInterface(method="PWPB95", basis="def2-TZVPP")
    inp_file = orca._generate_input(sample_xyz, tmp_path)

    content = inp_file.read_text()
    assert "def2-TZVPP/C" in content


def test_generate_input_multiple_double_hybrids(sample_xyz, tmp_path):
    """测试多种双杂化泛函"""
    double_hybrids = [
        ("PWPB95", "def2-TZVPP"),
        ("DSD-PBEP86", "def2-TZVPP"),
        ("B2PLYP", "def2-TZVPP"),
        ("DSD-BLYP", "def2-SVP")
    ]

    for method, basis in double_hybrids:
        orca = ORCAInterface(method=method, basis=basis)
        inp_file = orca._generate_input(sample_xyz, tmp_path)

        content = inp_file.read_text()
        # 验证 /C 辅助基组已添加
        c_basis = f"{basis}/C"
        assert c_basis in content
        assert method in content


def test_generate_input_non_double_hybrid(sample_xyz, tmp_path):
    """测试非双杂化泛函不添加 /C 基组"""
    # B3LYP 不是双杂化泛函
    orca = ORCAInterface(method="B3LYP", basis="def2-TZVPP")
    inp_file = orca._generate_input(sample_xyz, tmp_path)

    content = inp_file.read_text()
    # 不应该包含 /C 基组
    assert "def2-TZVPP/C" not in content
    assert "B3LYP" in content


# ==================== Session #2 测试用例 ====================

@pytest.fixture
def orca_output_success(tmp_path):
    """成功的 ORCA 输出文件"""
    out_file = tmp_path / "orca_success.out"
    out_file.write_text("""
-------------------
ORCA CALCULATION DONE
-------------------

Total Energy        :     -123.45678901
FINAL SINGLE POINT ENERGY      -123.45678901

ORCA TERMINATED NORMALLY
""")
    return out_file


@pytest.fixture
def orca_output_failed(tmp_path):
    """失败的 ORCA 输出文件"""
    out_file = tmp_path / "orca_failed.out"
    out_file.write_text("""
ORCA calculation encountered an error

ERROR: SCF failed to converge

ORCA TERMINATED WITH ERRORS
""")
    return out_file


@pytest.fixture
def orca_output_no_energy(tmp_path):
    """缺少能量信息的 ORCA 输出文件"""
    out_file = tmp_path / "orca_no_energy.out"
    out_file.write_text("""
ORCA calculation completed

ORCA TERMINATED NORMALLY
""")
    return out_file


def test_parse_output_energy(orca_output_success):
    """测试成功解析输出能量"""
    orca = ORCAInterface()
    result = orca._parse_output(orca_output_success)

    assert result.converged is True
    assert result.energy == -123.45678901
    assert result.error_message is None


def test_parse_output_failed(orca_output_failed):
    """测试解析失败输出"""
    orca = ORCAInterface()
    result = orca._parse_output(orca_output_failed)

    assert result.converged is False
    assert result.energy == 0.0
    assert result.error_message is not None
    assert "未正常终止" in result.error_message


def test_parse_output_no_energy(orca_output_no_energy):
    """测试缺少能量信息"""
    orca = ORCAInterface()
    result = orca._parse_output(orca_output_no_energy)

    assert result.converged is False
    assert result.energy == 0.0
    assert result.error_message is not None
    assert "无法找到能量信息" in result.error_message


# ==================== Session #3 测试用例 ====================

def test_find_orca_binary_with_path(tmp_path):
    """测试使用提供的路径查找 ORCA"""
    # 创建一个假的 ORCA 可执行文件
    fake_orca = tmp_path / "orca"
    fake_orca.write_text("#!/bin/bash\necho 'fake orca'")

    orca = ORCAInterface(orca_binary_path=str(fake_orca))

    assert orca.orca_binary is not None
    assert orca.orca_binary == fake_orca


def test_find_orca_binary_from_env():
    """测试从环境变量查找 ORCA"""
    test_path = '/tmp/test_orca_binary'
    with patch.dict('os.environ', {'ORCA_PATH': test_path}):
        with patch('shutil.which', return_value=test_path):
            with patch('pathlib.Path.exists', return_value=True):
                with patch('pathlib.Path.is_file', return_value=True):
                    orca = ORCAInterface()
                    assert orca.orca_binary is not None
                    assert str(orca.orca_binary) == test_path


def test_find_orca_binary_not_found():
    """测试未找到 ORCA 的情况"""
    # Mock 所有查找方式都失败
    with patch('shutil.which', return_value=None):
        with patch.dict('os.environ', {}, clear=True):
            orca = ORCAInterface()
            assert orca.orca_binary is None


@pytest.mark.skipif(
    True,  # 默认跳过，除非有 ORCA 环境
    reason="需要真实的 ORCA 环境"
)
def test_find_orca_binary_real():
    """测试查找真实的 ORCA (仅在 ORCA 可用时运行)"""
    orca = ORCAInterface()
    if orca.orca_binary is None:
        pytest.skip("ORCA 未找到")
    else:
        assert orca.orca_binary.exists()
        assert orca.orca_binary.is_file()


@patch('subprocess.Popen')
def test_run_orca_mock(mock_popen, sample_xyz, tmp_path):
    """测试运行 ORCA (Mock 版本)"""
    # 创建假的输入文件
    inp_file = tmp_path / "test.inp"
    inp_file.write_text("! M062X def2-SVP\n* xyz 0 1\nC 0 0 0\n*")

    # Mock 进程
    mock_process = MagicMock()
    mock_process.communicate.return_value = (None, None)
    mock_process.returncode = 0
    mock_popen.return_value = mock_process

    # 创建 ORCAInterface (设置假的可执行文件路径)
    fake_orca_path = tmp_path / "orca"
    fake_orca_path.write_text("fake orca")
    mpi_bin = tmp_path / "openmpi" / "bin"
    mpi_bin.mkdir(parents=True)
    (mpi_bin / "mpirun").write_text("#!/bin/sh\n")

    orca = ORCAInterface(orca_binary_path=str(fake_orca_path))

    # 运行 ORCA
    out_file = orca._run_orca(inp_file, tmp_path)

    # 验证
    assert out_file is not None
    assert out_file.suffix == '.out'
    mock_popen.assert_called_once()
    popen_kwargs = mock_popen.call_args.kwargs
    assert "env" in popen_kwargs
    assert isinstance(popen_kwargs["env"], dict)
    path_parts = popen_kwargs["env"]["PATH"].split(":")
    assert path_parts[0] == str(tmp_path)
    assert path_parts[1] == str(mpi_bin)
    assert "/opt/software/openmpi/bin" not in path_parts


def test_build_orca_runtime_env_allows_root_execution(tmp_path, monkeypatch):
    fake_orca_path = tmp_path / "orca"
    fake_orca_path.write_text("fake orca")
    mpi_bin = tmp_path / "openmpi" / "bin"
    mpi_bin.mkdir(parents=True)
    (mpi_bin / "mpirun").write_text("#!/bin/sh\n")
    monkeypatch.setenv("LD_LIBRARY_PATH", str(tmp_path / "openmpi" / "lib"))

    orca = ORCAInterface(orca_binary_path=str(fake_orca_path))
    env = orca._build_orca_runtime_env()

    assert env["ORCA_PATH"] == str(tmp_path)
    assert env["ORCA_DIR"] == str(tmp_path)
    assert env["OMPI_ALLOW_RUN_AS_ROOT"] == "1"
    assert env["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] == "1"


def test_generate_input_omits_pal_for_single_core(sample_xyz, tmp_path):
    orca = ORCAInterface(nprocs=1, orca_binary_path="/fake/orca")

    inp_file = orca._generate_input(sample_xyz, tmp_path)

    content = inp_file.read_text()
    assert "%pal nprocs" not in content


@patch('subprocess.Popen')
def test_run_orca_timeout(mock_popen, sample_xyz, tmp_path):
    """测试 ORCA 超时处理"""
    from subprocess import TimeoutExpired

    inp_file = tmp_path / "test.inp"
    inp_file.write_text("! M062X def2-SVP\n* xyz 0 1\nC 0 0 0\n*")

    # Mock 进程超时
    mock_process = MagicMock()
    mock_process.communicate.side_effect = TimeoutExpired('orca', 10)
    mock_process.kill = MagicMock()
    mock_popen.return_value = mock_process

    fake_orca_path = tmp_path / "orca"
    fake_orca_path.write_text("fake orca")

    orca = ORCAInterface(orca_binary_path=str(fake_orca_path))

    with pytest.raises(RuntimeError, match="超时"):
        orca._run_orca(inp_file, tmp_path, timeout=10)


def test_run_orca_no_binary(sample_xyz, tmp_path):
    """测试未找到 ORCA 二进制文件的情况"""
    inp_file = tmp_path / "test.inp"
    inp_file.write_text("! M062X def2-SVP\n* xyz 0 1\nC 0 0 0\n*")

    # 创建没有 ORCA 的接口
    orca = ORCAInterface()
    orca.orca_binary = None  # 强制设置为 None

    # 运行 ORCA，应该抛出 RuntimeError
    with pytest.raises(RuntimeError, match="ORCA 二进制文件未找到"):
        orca._run_orca(inp_file, tmp_path)


# ==================== Session #4 测试用例 ====================

@patch('rph_core.utils.orca_interface.ORCAInterface._run_orca')
def test_single_point_mock(mock_run_orca, sample_xyz, tmp_path):
    """测试端到端单点能计算 (Mock 版本)"""

    # Mock _run_orca 返回一个成功的输出文件
    out_file = tmp_path / "test.out"
    out_file.write_text("""
ORCA CALCULATION DONE
FINAL SINGLE POINT ENERGY      -234.56789012
ORCA TERMINATED NORMALLY
""")
    mock_run_orca.return_value = out_file

    # 创建假的 ORCA 可执行文件
    fake_orca = tmp_path / "orca"
    fake_orca.write_text("fake orca")

    orca = ORCAInterface(orca_binary_path=str(fake_orca))

    # 运行单点能计算
    result = orca.single_point(sample_xyz, tmp_path)

    # 验证结果
    assert result.converged is True
    assert result.energy == -234.56789012
    assert result.error_message is None

    # 验证 _run_orca 被调用
    mock_run_orca.assert_called_once()

    inp_files = list(tmp_path.glob("*.inp"))
    assert len(inp_files) >= 1, "Should generate at least one .inp file"


@patch('rph_core.utils.orca_interface.ORCAInterface._run_orca')
def test_single_point_mock_failed(mock_run_orca, sample_xyz, tmp_path):
    """测试端到端单点能计算失败 (Mock 版本)"""

    # Mock _run_orca 返回一个失败的输出文件
    out_file = tmp_path / "test.out"
    out_file.write_text("""
ORCA CALCULATION FAILED

ORCA TERMINATED WITH ERRORS
""")
    mock_run_orca.return_value = out_file

    fake_orca = tmp_path / "orca"
    fake_orca.write_text("fake orca")

    orca = ORCAInterface(orca_binary_path=str(fake_orca))

    # 运行单点能计算
    result = orca.single_point(sample_xyz, tmp_path)

    # 验证失败处理
    assert result.converged is False
    assert result.energy == 0.0
    assert result.error_message is not None
    assert "未正常终止" in result.error_message


@patch('rph_core.utils.orca_interface.ORCAInterface._run_orca')
def test_single_point_mock_exception(mock_run_orca, sample_xyz, tmp_path):
    """测试单点能计算异常处理 (Mock 版本)"""

    # Mock _run_orca 抛出异常
    mock_run_orca.side_effect = RuntimeError("ORCA 运行失败")

    fake_orca = tmp_path / "orca"
    fake_orca.write_text("fake orca")

    orca = ORCAInterface(orca_binary_path=str(fake_orca))

    # 运行单点能计算
    result = orca.single_point(sample_xyz, tmp_path)

    # 验证异常被捕获
    assert result.converged is False
    assert result.energy == 0.0
    assert result.error_message is not None
    assert "ORCA 运行失败" in result.error_message


# ==================== Session #7 测试用例 ====================

@pytest.mark.skipif(
    True,  # 默认跳过，除非有真实 ORCA 环境
    reason="需要真实的 ORCA 环境 - 设置 ORCA_PATH 环境变量以启用此测试"
)
def test_real_sp_orca_integration(sample_xyz, tmp_path):
    """
    测试真实 ORCA 单点能计算 (Session #7)

    注意: 此测试需要真实的 ORCA 环境才能运行
    要启用此测试:
    1. 安装 ORCA
    2. 设置环境变量: export ORCA_PATH=/path/to/orca
    3. 修改 @pytest.mark.skipif 条件为: os.environ.get('ORCA_PATH') is None
    """
    import os

    # 检查 ORCA 是否可用
    orca = ORCAInterface()
    if orca.orca_binary is None:
        pytest.skip("ORCA 二进制文件未找到")

    # 运行单点能计算
    result = orca.single_point(sample_xyz, tmp_path, timeout=60)

    # 验证结果
    assert result.converged is True, f"ORCA 计算未收敛: {result.error_message}"
    assert isinstance(result.energy, float), "能量应为浮点数"
    assert result.energy < 0, "能量应为负值 (Hartree)"

    # 验证输出文件存在
    out_file = tmp_path / "test.out"
    assert out_file.exists(), "ORCA 输出文件应存在"


@pytest.mark.skipif(
    True,  # 需要真实 ORCA 环境
    reason="需要真实的 ORCA 环境"
)
def test_real_sp_timeout_handling(sample_xyz, tmp_path):
    """
    测试真实 ORCA 超时处理 (Session #7)

    验证 ORCA 计算超时时能正确抛出 TimeoutError
    """
    import os

    orca = ORCAInterface()
    if orca.orca_binary is None:
        pytest.skip("ORCA 二进制文件未找到")

    # 设置非常短的超时时间 (1毫秒)
    with pytest.raises(TimeoutError, match="ORCA 计算超时"):
        orca.single_point(sample_xyz, tmp_path, timeout=1)


@pytest.mark.skipif(
    True,  # 需要真实 ORCA 环境
    reason="需要真实的 ORCA 环境"
)
def test_real_sp_invalid_xyz(tmp_path):
    """
    测试真实 ORCA 处理无效 XYZ (Session #7)

    验证 ORCA 能正确处理无效的输入文件
    """
    import os

    # 创建无效的 XYZ 文件
    invalid_xyz = tmp_path / "invalid.xyz"
    invalid_xyz.write_text("0\n\n")  # 空分子

    orca = ORCAInterface()
    if orca.orca_binary is None:
        pytest.skip("ORCA 二进制文件未找到")

    # 运行计算，应该失败
    result = orca.single_point(invalid_xyz, tmp_path, timeout=60)

    # 验证失败处理
    assert result.converged is False, "无效 XYZ 应导致计算失败"
    assert result.error_message is not None, "应有错误消息"


# ==================== Session #7 说明 ====================

"""
Session #7: 真实 ORCA 集成测试
================================

状态: ⏸️ 待执行 (需要真实 ORCA 环境)
优先级: P1 (可选)

完成内容:
- [x] 创建真实 ORCA 测试框架
- [x] 添加 test_real_sp_orca_integration() - 完整单点能计算
- [x] 添加 test_real_sp_timeout_handling() - 超时处理
- [x] 添加 test_real_sp_invalid_xyz() - 无效输入处理

验收标准:
✅ 测试框架已创建，包含 3 个测试用例
✅ 所有测试默认跳过 (需要 ORCA 环境)
✅ 测试包含详细的启用说明

如何启用:
1. 安装 ORCA 量子化学软件
2. 设置环境变量: export ORCA_PATH=/path/to/orca
3. 修改 @pytest.mark.skipif 条件
4. 运行: pytest tests/test_orca_interface.py::test_real_sp -v

预期结果 (有 ORCA 时):
- test_real_sp_orca_integration: 成功返回负能量值
- test_real_sp_timeout_handling: 正确抛出 TimeoutError
- test_real_sp_invalid_xyz: 返回 converged=False

阻塞原因:
- 需要 ORCA 软件 (可能需要学术许可)
- 计算时间较长 (不适合 CI/CD)

替代方案:
- Mock 测试 (Sessions #1-#6 已完成)
- 验证脚本 (verify_orca.py, verify_session2.py 已创建)
"""
