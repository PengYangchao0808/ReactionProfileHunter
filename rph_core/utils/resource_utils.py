"""
Resource Utilities Module
=======================
资源管理工具：内存换算、路径解析、环境变量处理

Author: QC Descriptors Team
Date: 2026-01-13
Purpose: 替代旧版 Bash 脚本中的资源管理逻辑
"""

import re
import os
import shutil
import logging
from pathlib import Path
from typing import Optional, Dict, Tuple

logger = logging.getLogger(__name__)


def get_project_root() -> Path:
    """
    获取项目根目录的绝对路径 (本文件 rph_core/utils/resource_utils.py 的上级两级)

    Returns:
        项目根目录的绝对路径
    """
    current_file = Path(__file__).resolve()
    project_root = current_file.parent.parent.parent
    return project_root


# ========== 内存换算函数 ==========
def mem_to_mb(mem_str: str) -> int:
    """
    将内存字符串转换为 MB

    支持格式:
    - "32GB", "16GB", "32G" -> 32 * 1024 MB
    - "4096MB", "2048MB", "4096M" -> 直接返回 MB 值
    - 纯数字 (如 "16") -> 视为 GB

    Args:
        mem_str: 内存字符串

    Returns:
        内存大小（MB）
    """
    if not mem_str:
        logger.warning("内存字符串为空，返回默认值 4000 MB")
        return 4000

    m = mem_str.strip().upper()
    match = re.match(r'^(\d+)\s*(GB?|MB?)?$', m)

    if match:
        val = int(match.group(1))
        unit = match.group(2) or 'GB'

        if unit.startswith('G'):
            return val * 1024
        else:
            return val

    logger.warning(f"无法解析内存字符串: {mem_str}，返回默认值 4000 MB")
    return 4000


def calc_orca_maxcore(mem: str, nproc: int, safety_factor: float = 0.65) -> int:
    """
    计算 ORCA 的 maxcore 参数（MB per core）

    旧脚本逻辑（Config.sh 第90-96行）:
    ```bash
    if [[ "$MEM" =~ ^([0-9]+)[[:space:]]*GB$ ]]; then
      ORCA_MAXCORE_MB="$(( ${BASH_REMATCH[1]} * 1024 * 2 ))"
    fi
    ```

    Args:
        mem: 内存字符串（如 "32GB"）
        nproc: 并行核数
        safety_factor: 安全系数（避免内存溢出）

    Returns:
        maxcore 值（MB per core）
    """
    total_mb = mem_to_mb(mem)
    maxcore = int(total_mb * safety_factor / nproc)

    logger.debug(f"ORCA maxcore 计算: {mem} / {nproc} cores * {safety_factor} = {maxcore} MB/core")
    return maxcore


# ========== 路径解析函数 ==========
def find_executable(
    program_name: str,
    config_path: Optional[str] = None,
    env_vars: Optional[list] = None,
    allow_path_search: bool = True,
    known_dirs: Optional[list] = None,
    binary_names: Optional[list] = None,
) -> Optional[Path]:
    """
    查找可执行文件路径（多级回退策略）

    查找顺序:
    1. 配置文件中的绝对路径
    2. 环境变量
    3. 系统 PATH
    4. 已知安装目录扫描

    Args:
        program_name: 程序名称（如 "orca"，用于 PATH 搜索）
        config_path: 配置文件中的路径
        env_vars: 要检查的环境变量列表
        allow_path_search: 是否允许在系统 PATH 中查找
        known_dirs: 已知安装目录列表
        binary_names: 真实 binary 名称列表（如 ["g16"]），用于 known_dirs 搜索

    Returns:
        可执行文件的 Path 对象，如果找不到返回 None
    """
    # 1. 配置文件路径
    if config_path:
        path = Path(config_path)
        if path.exists() and path.is_file() and os.access(path, os.X_OK):
            logger.debug(f"Found: {program_name} (config): {path}")
            return path
        else:
            logger.warning(f"Config path invalid: {config_path}")

    # 2. 环境变量
    if env_vars:
        for env_var in env_vars:
            env_path = os.environ.get(env_var)
            if env_path:
                path = Path(env_path)
                if path.exists() and path.is_file() and os.access(path, os.X_OK):
                    logger.debug(f"Found: {program_name} (env {env_var}): {path}")
                    return path
                if path.is_dir():
                    names_to_try = binary_names or [program_name]
                    for name in names_to_try:
                        for candidate in (path / name, path / 'bin' / name):
                            if candidate.is_file() and os.access(str(candidate), os.X_OK):
                                logger.debug(f"Found: {program_name} (env {env_var} dir): {candidate}")
                                return candidate
                logger.warning(f"Env {env_var}={env_path} points to non-existent or invalid file/dir")

    # 3. 系统 PATH
    if allow_path_search:
        which_result = shutil.which(program_name)
        if which_result:
            path = Path(which_result)
            logger.debug(f"Found: {program_name} (PATH): {path}")
            return path

    # 4. 已知安装目录扫描
    if known_dirs:
        names_to_try = binary_names or [program_name]
        import glob as _glob
        for pattern in known_dirs:
            for base in _glob.glob(pattern):
                for name in names_to_try:
                    # Try: base/name
                    candidate = Path(base) / name
                    if candidate.is_file() and os.access(str(candidate), os.X_OK):
                        logger.debug(f"Found: {program_name} (known_dirs): {candidate}")
                        return candidate
                    # Try: base/bin/name
                    candidate = Path(base) / 'bin' / name
                    if candidate.is_file() and os.access(str(candidate), os.X_OK):
                        logger.debug(f"Found: {program_name} (known_dirs/bin): {candidate}")
                        return candidate

    logger.warning(f"Not found: {program_name}. Searched: config={config_path}, env={env_vars}, PATH, known_dirs={known_dirs}")
    return None


def setup_ld_library_path(ld_library_paths: list) -> Optional[str]:
    """
    设置 LD_LIBRARY_PATH 环境变量（ORCA MPI 支持）

    幂等：已存在的路径条目不会被重复追加。

    Args:
        ld_library_paths: 要添加的路径列表

    Returns:
        更新后的完整 LD_LIBRARY_PATH 字符串；无变化时返回 None。
    """
    if not ld_library_paths:
        return None

    current_path = os.environ.get('LD_LIBRARY_PATH', '')
    existing_entries = set(current_path.split(':')) if current_path else set()
    new_paths = []

    for path_str in ld_library_paths:
        individual_paths = path_str.split(':') if ':' in path_str else [path_str]
        for single_path in individual_paths:
            single_path = single_path.strip()
            if not single_path:
                continue
            path = Path(single_path)
            resolved = str(path)
            if path.exists() and resolved not in existing_entries:
                new_paths.append(resolved)
                existing_entries.add(resolved)
            elif not path.exists():
                logger.warning(f"库路径不存在: {single_path}")

    if not new_paths:
        return None

    if current_path:
        new_path = ":".join(new_paths) + ":" + current_path
    else:
        new_path = ":".join(new_paths)

    os.environ['LD_LIBRARY_PATH'] = new_path
    logger.info(f"已更新 LD_LIBRARY_PATH: {new_path}")
    return new_path


# ========== 配置派生函数 ==========
def resolve_executable_config(
    config: dict,
    program_key: str,
    env_vars: Optional[list] = None,
    allow_path_search: bool = True,
    known_dirs: Optional[list] = None,
    fail_on_invalid: bool = True,
    binary_names: Optional[list] = None,
) -> Dict:
    exe_cfg = (config.get('executables', {}) or {}).get(program_key, {})
    if not isinstance(exe_cfg, dict):
        exe_cfg = {}
    res_cfg = config.get('resources', {}) or {}

    # 检查显式路径是否存在（或是否应 fail-fast）
    config_path = exe_cfg.get('path')
    if config_path:
        p = Path(config_path)
        if not p.exists() or not p.is_file():
            if fail_on_invalid:
                raise RuntimeError(
                    f"Config executables.{program_key}.path='{config_path}' "
                    f"not found and fail_on_invalid_explicit=True. "
                    f"Set the correct path or enable auto-discovery."
                )
            else:
                logger.warning(
                    f"Config path '{config_path}' invalid; falling back to auto-discovery"
                )

    exe_path = find_executable(
        program_key,
        config_path=config_path,
        env_vars=env_vars,
        allow_path_search=allow_path_search,
        known_dirs=known_dirs,
        binary_names=binary_names,
    )

    source = "not_found"
    if exe_path:
        if config_path and exe_path == Path(config_path):
            source = "config"
        elif env_vars:
            for env_var in env_vars:
                env_path = os.environ.get(env_var)
                if env_path and exe_path == Path(env_path):
                    source = "env"
                    break
            if source == "not_found":
                which_result = shutil.which(program_key)
                if which_result and exe_path == Path(which_result):
                    source = "path"
                else:
                    source = "known_dirs"
        else:
            source = "path"

    ld_path_str = exe_cfg.get('ld_library_path')

    if ld_path_str:
        setup_ld_library_path([ld_path_str])

    return {
        'path': exe_path,
        'ld_library_path': ld_path_str,
        'found': exe_path is not None,
        'source': source
    }


# ═══════════════════════════════════════════════════
# 可执行文件发现常量表
# ═══════════════════════════════════════════════════

_BINARY_NAMES = {
    'gaussian': 'g16',
    'orca': 'orca',
    'xtb': 'xtb',
    'crest': 'crest',
    'isostat': 'isostat',
    'shermo': 'Shermo',
}

_KNOWN_DIRS = {
    'gaussian': ['/opt/software/gaussian', '/opt/g16', '/opt/gaussian', '/usr/local/g16'],
    'orca': ['/opt/software/orca', '/opt/orca', '/usr/local/orca'],
    'xtb': ['/opt/software/xtb', '/opt/xtb', '/usr/local/bin'],
    'crest': ['/opt/software/crest', '/opt/crest', '/usr/local/bin'],
    'isostat': ['/opt/software/molclus', '/opt/molclus', '/usr/local/bin', '/usr/bin'],
    'shermo': ['/opt/software/shermo', '/opt/shermo', '/usr/local/bin'],
}

_ENV_VARS = {
    'gaussian': ['GAUSS_PATH', 'GAUSSIAN_PATH'],
    'orca': ['ORCA_PATH', 'ORCA_BIN'],
    'xtb': ['XTB_PATH', 'XTB_BIN'],
    'crest': ['CREST_PATH', 'CREST_BIN'],
    'isostat': ['ISOSTAT_PATH', 'MOLCLUS_HOME'],
    'shermo': ['SHERMO_PATH'],
}

_MPI_KNOWN_DIRS = ['/opt/openmpi', '/opt/software/openmpi', '/usr/lib/openmpi', '/usr/lib64/openmpi']
_MPI_ENV_VARS = ['MPI_BIN_DIR', 'MPI_HOME']


def resolve_all_executables(config: dict) -> Dict[str, Dict]:
    """一次解析全部可执行文件，返回 source 报告。

    在 orchestrator 启动时调用一次，位于 resolve_resources() 之后。
    根据运行时上下文判断各程序是否为必需项。
    """
    programs = ['gaussian', 'orca', 'xtb', 'crest', 'isostat', 'shermo']
    discovery_cfg = (config.get('executables', {}) or {}).get('discovery', {})
    if not isinstance(discovery_cfg, dict):
        discovery_cfg = {}
    fail_on_invalid = bool(discovery_cfg.get('fail_on_invalid_explicit', True))
    allow_discovery = bool(discovery_cfg.get('enabled', True))
    yaml_known_dirs = discovery_cfg.get('known_dirs', {})
    if not isinstance(yaml_known_dirs, dict):
        yaml_known_dirs = {}

    results = {}
    for prog in programs:
        if allow_discovery:
            known = yaml_known_dirs.get(prog) or _KNOWN_DIRS.get(prog)
        else:
            known = None
        results[prog] = resolve_executable_config(
            config, prog,
            env_vars=_ENV_VARS.get(prog),
            known_dirs=known,
            binary_names=[_BINARY_NAMES.get(prog, prog)],
            fail_on_invalid=fail_on_invalid,
        )

    # MPI（ORCA 辅助程序）
    results['mpirun'] = _resolve_mpi(config, allow_discovery)

    return results


def _resolve_mpi(config: dict, allow_discovery: bool) -> Dict:
    """解析 MPI (mpirun) 路径。"""
    import glob as _glob

    exe_root = config.setdefault('executables', {})
    if not isinstance(exe_root, dict):
        exe_root = {}
        config['executables'] = exe_root
    orca_cfg = exe_root.setdefault('orca', {})
    if not isinstance(orca_cfg, dict):
        orca_cfg = {}
        exe_root['orca'] = orca_cfg
    mpi_bin = orca_cfg.get('mpi_bin_dir')
    if mpi_bin and Path(mpi_bin).is_dir() and (Path(mpi_bin) / 'mpirun').exists():
        return {'path': Path(mpi_bin) / 'mpirun', 'found': True, 'source': 'config'}

    # 环境变量
    for var in _MPI_ENV_VARS:
        val = os.environ.get(var)
        if val:
            candidate = Path(val) / 'mpirun'
            if candidate.is_file():
                orca_cfg['mpi_bin_dir'] = str(Path(val))
                orca_cfg['mpi_lib_dir'] = str(Path(val).parent / 'lib')
                return {'path': candidate, 'found': True, 'source': 'env'}

    # PATH
    which_mpi = shutil.which('mpirun')
    if which_mpi:
        mpi_dir = str(Path(which_mpi).parent)
        orca_cfg.setdefault('mpi_bin_dir', mpi_dir)
        orca_cfg.setdefault('mpi_lib_dir', str(Path(mpi_dir).parent / 'lib'))
        return {'path': Path(which_mpi), 'found': True, 'source': 'path'}

    # known_dirs
    if allow_discovery:
        for pattern in _MPI_KNOWN_DIRS:
            for base in _glob.glob(pattern):
                mpirun = Path(base) / 'bin' / 'mpirun'
                if mpirun.is_file():
                    orca_cfg.setdefault('mpi_bin_dir', str(mpirun.parent))
                    orca_cfg.setdefault('mpi_lib_dir', str(Path(base) / 'lib'))
                    return {'path': mpirun, 'found': True, 'source': 'known_dirs'}

    return {'path': None, 'found': False, 'source': 'not_found'}


def _is_executable_required(config: dict, program_key: str) -> bool:
    """根据配置判断该可执行文件是否必需。"""
    theory = config.get('theory', {})
    step1 = config.get('step1', {})

    REQUIRED_MAP = {
        'gaussian': (
            theory.get('optimization', {}).get('engine') == 'gaussian'
            or theory.get('single_point', {}).get('engine') == 'gaussian'
        ),
        'orca': (
            theory.get('optimization', {}).get('engine') == 'orca'
            or theory.get('single_point', {}).get('engine') == 'orca'
        ),
        'xtb': True,
        'crest': not step1.get('conformer_search', {}).get('two_stage_enabled', False),
        'isostat': True,
        'shermo': True,
        'mpirun': (
            theory.get('single_point', {}).get('engine') == 'orca'
            and (theory.get('single_point', {}).get('nproc') or 1) > 1
        ),
    }
    return REQUIRED_MAP.get(program_key, False)


def declare_resolved_executables(config: dict, results: Dict[str, Dict]) -> None:
    """将 resolve_all_executables() 的结果回写到 config['executables'][prog]['path']。

    仅当 found=True 且当前 config 中的 path 缺失或无效时才回写，
    避免覆盖用户显式指定的有效路径。
    """
    exe_root = config.setdefault('executables', {})
    if not isinstance(exe_root, dict):
        exe_root = {}
        config['executables'] = exe_root

    for prog, info in results.items():
        if prog == 'mpirun':
            continue
        if not info.get('found'):
            continue
        path = info.get('path')
        if path is None:
            continue

        prog_cfg = exe_root.setdefault(prog, {})
        if not isinstance(prog_cfg, dict):
            prog_cfg = {}
            exe_root[prog] = prog_cfg

        current_path = prog_cfg.get('path')
        if current_path and Path(current_path).exists():
            continue

        prog_cfg['path'] = str(path)
        source = info.get('source', 'unknown')
        if source != 'config':
            logger.debug(f"[declare] executables.{prog}.path = {path} (source={source})")


def get_gaussian_command(config: dict) -> Tuple[str, bool]:
    """统一解析 Gaussian 命令路径。

    All Gaussian execution paths should use this function,
    避免硬编码 ``"g16"`` 或忽略 ``use_wrapper`` 开关。

    Returns:
        (command, use_wrapper) — command 是可调用的命令字符串；
        use_wrapper 表示是否通过 wrapper 脚本调用。
    """
    exe_cfg = (config.get('executables', {}) or {}).get('gaussian', {})
    if not isinstance(exe_cfg, dict):
        exe_cfg = {}
    use_wrapper = bool(exe_cfg.get('use_wrapper', True))

    if use_wrapper:
        wrapper_path = exe_cfg.get('wrapper_path', './scripts/run_g16_worker.sh')
        wrapper_full = Path(wrapper_path)
        if not wrapper_full.is_absolute():
            wrapper_full = (get_project_root() / wrapper_path).resolve()
        if not wrapper_full.exists():
            raise RuntimeError(
                f"Gaussian wrapper not found: {wrapper_full}. "
                f"Set executables.gaussian.use_wrapper=false or fix wrapper_path."
            )
        return str(wrapper_full), True

    resolved = resolve_executable_config(
        config, 'gaussian', env_vars=['GAUSS_PATH', 'GAUSSIAN_PATH']
    )
    path = resolved.get('path')
    if path:
        return str(path), False
    logger.warning("Gaussian executable not resolved; falling back to literal 'g16'")
    return 'g16', False


# ========== 验证函数 ==========
def validate_executable(path: Path, program_name: str = "") -> bool:
    """
    验证可执行文件是否可用

    Args:
        path: 可执行文件路径
        program_name: 程序名称（用于错误提示）

    Returns:
        是否可用
    """
    if not path:
        logger.error(f"{program_name}: 路径为空")
        return False

    if not path.exists():
        logger.error(f"{program_name}: 文件不存在: {path}")
        return False

    if not path.is_file():
        logger.error(f"{program_name}: 不是文件: {path}")
        return False

    if not os.access(path, os.X_OK):
        logger.error(f"{program_name}: 文件不可执行: {path}")
        return False

    return True


def _nested_get(config: dict, key_path: str):
    """从嵌套 dict 中按点号分隔的 key_path 取值。"""
    parts = key_path.split('.')
    d = config
    for part in parts:
        if not isinstance(d, dict):
            return None
        d = d.get(part)
        if d is None:
            return None
    return d


def _nested_set(config: dict, key_path: str, value):
    """在嵌套 dict 中按点号分隔的 key_path 设置值。"""
    parts = key_path.split('.')
    d = config
    for part in parts[:-1]:
        d = d.setdefault(part, {})
    d[parts[-1]] = value


def _auto_detect_mem_gb() -> int:
    """自动检测系统物理内存 (GB)，失败返回 32。"""
    try:
        import psutil
        return max(8, int(psutil.virtual_memory().total / (1024**3) * 0.75))
    except ImportError:
        pass
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemTotal:'):
                    total_kb = int(line.split()[1])
                    return max(8, int(total_kb / (1024 * 1024) * 0.75))
    except (IOError, ValueError):
        pass
    return 32


def resolve_resources(config: dict) -> dict:
    """
    统一资源配置解析函数。

    在 pipeline 启动时调用一次，执行以下操作：
    1. 若 resources.* 为 null → 自动检测系统 CPU/内存
    2. 各子区段（theory.*, step*.*, intra_reaction_parallel.*）若为 null → 继承
    3. 并行 lane 做预算分配：per_job = total_cores / workers
    4. 显式数值保持不动（override 语义）
    """
    resources = config.setdefault('resources', {})

    # ── Step 1: 自动检测系统资源 ──
    if not resources.get('nproc'):
        detected = os.cpu_count() or 1
        resources['nproc'] = max(1, detected)
        logger.info(f"Auto-detected CPU: {detected} cores → nproc={resources['nproc']}")

    if not resources.get('mem'):
        detected_gb = _auto_detect_mem_gb()
        resources['mem'] = f"{detected_gb}GB"
        logger.info(f"Auto-detected RAM: using {resources['mem']}")

    if not resources.get('orca_maxcore_safety'):
        resources['orca_maxcore_safety'] = 0.65

    nproc = resources['nproc']
    mem = resources['mem']

    # ── Step 2: 继承 resources.* 到 theory.* ──
    _inherit_if_null(config, 'theory.optimization.nproc', nproc)
    _inherit_if_null(config, 'theory.optimization.mem', mem)
    _inherit_if_null(config, 'theory.single_point.nproc', nproc)

    # ── Step 3: 继承到 step1 ──
    _inherit_if_null(config, 'step1.crest.threads', nproc)
    _inherit_if_null(config, 'step1.conformer_search.common.threads',
                     _nested_get(config, 'step1.crest.threads') or nproc)

    # ── Step 4: 继承到 step2 ──
    _inherit_if_null(config, 'step2.scan.nproc', nproc)
    _inherit_if_null(config, 'step2.xtb_settings.nproc', nproc)
    _inherit_if_null(config, 'step2.crest_rescue.threads', nproc)

    # ── Step 5: 并行预算分配 ──
    # Always resolve the optional parallel budget from the global resources.
    irp = config.get('intra_reaction_parallel', {})
    _resolve_irp_defaults(irp, resources, nproc, mem)

    return config


def _resolve_irp_defaults(irp: dict, resources: dict, nproc: int, mem: str):
    """填充 intra_reaction_parallel 区段的默认值。

    始终执行（不受 enabled 标志控制），确保所有 lane 有合理默认值。
    并行预算原则：多个 lane 共享 total_cores/total_mem，不得超分配。
    """
    from rph_core.utils.resource_utils import mem_to_mb

    total_cores = irp.get('total_cores') or nproc
    total_mem_str = irp.get('total_mem') or mem
    total_mem_mb = mem_to_mb(total_mem_str)
    safety = resources.get('orca_maxcore_safety', 0.65)

    # 确保子结构存在
    irp.setdefault('s1', {})
    irp.setdefault('s2', {})
    irp.setdefault('s3', {})
    s1 = irp['s1']
    s2 = irp['s2']
    s3 = irp['s3']

    # ── S1 分子 lane ──
    mol_workers = s1.get('max_molecule_workers', 1)
    conf_workers = s1.get('max_conformer_workers', 1)

    if s1.get('crest_cores') is None:
        s1['crest_cores'] = max(1, int(total_cores / max(mol_workers, 1)))
    if s1.get('opt_cores_per_job') is None:
        s1['opt_cores_per_job'] = max(1, int(total_cores / max(mol_workers, 1)))
    if s1.get('molecule_mem_per_job') is None:
        s1['molecule_mem_per_job'] = f"{max(8, int(total_mem_mb / max(mol_workers, 1) / 1024))}GB"
    if s1.get('opt_mem_per_job') is None:
        s1['opt_mem_per_job'] = s1['molecule_mem_per_job']

    if s1.get('sp_cores_per_job') is None:
        s1['sp_cores_per_job'] = max(1, int(total_cores / max(conf_workers, 1)))
    if s1.get('sp_maxcore_per_job') is None:
        s1['sp_maxcore_per_job'] = int(total_mem_mb * safety / max(conf_workers, 1) / s1['sp_cores_per_job'])

    # ── S2 ──
    if s2.get('nproc') is None:
        s2['nproc'] = total_cores

    # ── S3 ──
    parallel_it = s3.get('parallel_intermediate_ts', True)
    sp_workers = s3.get('sp_workers', 1)

    if s3.get('ts_opt_cores') is None:
        if parallel_it:
            s3['ts_opt_cores'] = max(1, int(total_cores * 3 / 4))
        else:
            s3['ts_opt_cores'] = total_cores

    if s3.get('intermediate_opt_cores') is None:
        if parallel_it:
            s3['intermediate_opt_cores'] = max(1, int(total_cores / 4))
        else:
            s3['intermediate_opt_cores'] = total_cores

    if s3.get('ts_opt_mem') is None:
        if parallel_it:
            ts_mem = max(8, int(total_mem_mb * 3 / 4 / 1024))
            s3['ts_opt_mem'] = f"{ts_mem}GB"
        else:
            s3['ts_opt_mem'] = total_mem_str

    if s3.get('intermediate_opt_mem') is None:
        if parallel_it:
            inter_mem = max(4, int(total_mem_mb / 4 / 1024))
            s3['intermediate_opt_mem'] = f"{inter_mem}GB"
        else:
            s3['intermediate_opt_mem'] = total_mem_str

    if s3.get('sp_cores_per_job') is None:
        s3['sp_cores_per_job'] = max(1, int(total_cores / max(sp_workers, 1)))
    if s3.get('sp_maxcore_per_job') is None:
        s3['sp_maxcore_per_job'] = int(
            total_mem_mb * safety / max(sp_workers, 1) / s3['sp_cores_per_job']
        ) if s3['sp_cores_per_job'] > 0 else 4000

    # ── 回写 total_cores/total_mem ──
    irp['total_cores'] = total_cores
    irp['total_mem'] = total_mem_str

    # ── 最终预算校验 ──
    _validate_irp_budget(irp, resources, total_cores, total_mem_mb)


def _validate_irp_budget(irp: dict, resources: dict, total_cores: int, total_mem_mb: int):
    """校验所有 lane 合并后不超总预算（覆盖自动推导和显式 override 的最终值）。"""
    logger = logging.getLogger(__name__)
    safety = resources.get('orca_maxcore_safety', 0.65)

    s3 = irp.get('s3', {})
    parallel_it = s3.get('parallel_intermediate_ts', True)

    if parallel_it:
        ts_cores = s3.get('ts_opt_cores', 0) or 0
        inter_cores = s3.get('intermediate_opt_cores', 0) or 0
        if ts_cores + inter_cores > total_cores:
            logger.warning(
                f"S3 budget: TS({ts_cores})+Inter({inter_cores})="
                f"{ts_cores+inter_cores} > total({total_cores})"
            )

    for lane_name in ('s1', 's3'):
        lane = irp.get(lane_name, {})
        sp_cores = lane.get('sp_cores_per_job', 0) or 0
        sp_maxcore = lane.get('sp_maxcore_per_job', 0) or 0
        sp_workers = lane.get('sp_workers', 1)
        if sp_cores > 0 and sp_maxcore > 0:
            total_sp_mem = sp_cores * sp_maxcore * sp_workers
            budget_mem = int(total_mem_mb * safety)
            if total_sp_mem > budget_mem:
                logger.warning(
                    f"{lane_name} SP: {total_sp_mem/1024:.1f}GB > "
                    f"budget {budget_mem/1024:.1f}GB"
                )


def _inherit_if_null(config: dict, key_path: str, default):
    """若 key_path 值为 None/null/0/""，设为 default。"""
    parts = key_path.split('.')
    d = config
    for part in parts[:-1]:
        d = d.setdefault(part, {})
    val = d.get(parts[-1])
    if val is None or val == 0 or (isinstance(val, str) and not val.strip()):
        d[parts[-1]] = default
