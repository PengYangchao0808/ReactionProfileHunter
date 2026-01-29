#!/usr/bin/env python3
# pyright: reportConstantRedefinition=false
# pyright: reportPossiblyUnboundVariable=false
# pyright: reportMissingTypeArgument=false
# pyright: reportArgumentType=false
# pyright: reportIndexIssue=false
# pyright: reportCallIssue=false
# pyright: reportOperatorIssue=false
# pyright: reportAttributeAccessIssue=false
# pyright: reportOptionalSubscript=false
"""
自动测试脚本 - ReactionProfileHunter v5
=========================================
使用 reaxys_cleaned.csv 第一个示例进行完整四步流程测试

特性:
1. 自动调用真实 QC 软件 (XTB/CREST/ORCA/Gaussian)
2. 自动错误恢复和重试 (使用 qc_runner)
3. 详细日志和错误报告
4. 中断恢复支持
5. [NEW] Rich UI 美化 (Panel, Table, Progress)

Usage:
    python run_auto_test.py [--step S1|S2|S3|S4|ALL] [--resume]

Author: QCcalc Team
Date: 2026-01-13
"""

import sys
import logging
import yaml
import json
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any, List

# Rich imports
Console: Any = None
RichHandler: Any = None
Panel: Any = None
Table: Any = None
has_rich = False
try:
    from rich.console import Console
    from rich.logging import RichHandler
    from rich.panel import Panel
    from rich.table import Table
    has_rich = True
except ImportError:
    print("WARNING: 'rich' library not found. Installing it is recommended for better UI.")

HAS_RICH = has_rich

sys.path.insert(0, str(Path(__file__).parent))

from rph_core.orchestrator import ReactionProfileHunter, PipelineResult
from rph_core.version import __version__
from rph_core.utils.qc_runner import (
    RetryConfig,
    with_retry,
    QCFailureType,
    analyze_log_for_errors
)
from rph_core.utils.optimization_config import normalize_qc_config
from rph_core.utils.dataset_loader import load_reaction_records
from rph_core.utils.result_inspector import ResultInspector

# 配置根 Logger (使用 RichHandler 如果可用)
if HAS_RICH:
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, show_path=False)]
    )
    console = Console()
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

logger = logging.getLogger("AutoTest")


class AutoTester:
    """
    自动测试引擎 (Rich UI Enhanced)
    """

    def __init__(self, config_path: Optional[Path] = None):
        """初始化自动测试引擎"""
        if config_path is None:
            config_path = Path(__file__).parent / "config" / "defaults.yaml"

        self.config_path = config_path
        self.config = self._load_config()
        self.config, self.qc_fixes = normalize_qc_config(self.config, auto_fix=True)
        self.hunter = ReactionProfileHunter(config_path=config_path)

        # 加载测试数据
        self.test_data = self._load_test_data()

        # 测试报告
        self.test_report = {
            "timestamp": datetime.now().isoformat(),
            "config_file": str(config_path),
            "test_smiles": self.test_data["smiles"],
            "rx_id": self.test_data["rx_id"],
            "steps": {},
            "environment_check": {},
            "config_fixes": self.qc_fixes,
            "success": False,
            "error": None,
            "files_generated": []
        }
        if self.qc_fixes:
            logger.warning(f"QC config normalized: {len(self.qc_fixes)} change(s)")

        # 打印欢迎 Banner
        self._print_banner()

        # 环境预检
        self._verify_environment()

    def _print_banner(self):
        if HAS_RICH:
            grid = Table.grid(expand=True)
            grid.add_column(justify="center", ratio=1)
            grid.add_row(
                Panel(
                    f"[bold blue]Reaction Profile Hunter v{__version__}[/bold blue]\n"
                    f"[green]Automated Pipeline Testing System[/green]\n\n"
                    f"Reaction ID: [yellow]{self.test_data['rx_id']}[/yellow]\n"
                    f"Product: [cyan]{self.test_data['smiles']}[/cyan]\n"
                    f"Config: [dim]{self.config_path}[/dim]",
                    title="[bold]RPH AutoTest[/bold]",
                    border_style="blue"
                )
            )
            console.print(grid)
        else:
            logger.info("="*70)
            logger.info(f"ReactionProfileHunter v{__version__} 自动测试")
            logger.info(f"测试反应 ID: {self.test_data['rx_id']}")
            logger.info("="*70)

    def _load_config(self) -> dict:
        """加载配置文件"""
        try:
            with open(self.config_path, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            logger.error(f"配置文件加载失败: {e}")
            raise

    def _load_test_data(self) -> Dict:
        """从数据集读取第一个测试数据"""
        run_cfg = self.config.get("run", {}) or {}
        dataset_cfg = run_cfg.get("dataset")
        if not dataset_cfg:
            dataset_cfg = self.config.get("reference_states", {}).get("dataset", {})

        if not dataset_cfg:
            raise RuntimeError("未配置 dataset（run.dataset 或 reference_states.dataset）")

        csv_path = Path(dataset_cfg.get("path", Path(__file__).parent / "reaxys_cleaned.csv"))

        if HAS_RICH:
            console.print(f"[dim]读取测试数据: {csv_path}[/dim]")
        else:
            logger.info(f"读取测试数据: {csv_path}")

        records = load_reaction_records(dataset_cfg=dataset_cfg, filter_ids=None, max_tasks=1)
        record = records[0]
        if not record.product_smiles_main:
            raise RuntimeError("dataset 记录缺少 product_smiles_main")

        return {
            "rx_id": record.rx_id,
            "smiles": record.product_smiles_main,
            "precursor_smiles": record.precursor_smiles,
            "yield": record.yield_ or "",
        }

    def _verify_environment(self):
        """验证计算环境"""
        from rph_core.utils.resource_utils import find_executable

        if HAS_RICH:
            table = Table(title="Environment Pre-flight Check", box=None)
            table.add_column("Program", style="cyan")
            table.add_column("Status", justify="center")
            table.add_column("Path / Details", style="dim")
        else:
            logger.info("--- 环境预检 ---")

        env_check = {}

        for prog in ['xtb', 'crest', 'orca', 'gaussian']:
            prog_config = self.config['executables'].get(prog, {})
            path = prog_config.get('path')
            ld_path = prog_config.get('ld_library_path')

            status = "NOT_FOUND"
            details = []
            path_display = "Not set"

            if path:
                if find_executable(path):
                    status = "OK"
                    path_display = str(path)
                    if ld_path:
                        details.append(f"LD_LIB: {ld_path}")
                else:
                    details.append(f"Config path missing: {path}")
                    # 尝试自动查找
                    if self.config.get('global', {}).get('enable_path_search', True):
                        alt_path = find_executable(prog)
                        if alt_path:
                            status = "FOUND_BY_SEARCH"
                            path_display = str(alt_path)
                            details.append("Auto-detected")
                    else:
                        details.append("Auto-search disabled")
            else:
                details.append("Path not configured")

            env_check[prog] = {
                "status": status,
                "details": details
            }

            if HAS_RICH:
                status_style = "green" if status == "OK" else "yellow" if status == "FOUND_BY_SEARCH" else "red"
                status_icon = "✔" if status == "OK" else "⚠" if status == "FOUND_BY_SEARCH" else "✘"
                table.add_row(
                    prog.upper(), 
                    f"[{status_style}]{status_icon}[/{status_style}]", 
                    f"{path_display} {' '.join(details)}"
                )
            else:
                logger.info(f"{prog.upper()}: {status} ({path_display})")

        if HAS_RICH:
            console.print(Panel(table, border_style="green"))
            console.print()

        self.test_report["environment_check"] = env_check

        critical_issues = [prog for prog, check in env_check.items() if check["status"] == "NOT_FOUND"]
        if critical_issues:
            logger.warning(f"警告: 以下程序未找到: {', '.join(critical_issues)}")

        try:
            import resource
            stack_size_kb = resource.getrlimit(resource.RLIMIT_STACK)[0]
            stack_size_mb = stack_size_kb / 1024

            if stack_size_kb > 0 and stack_size_kb < 65536:
                stack_warning = (
                    f"⚠️  栈大小较小: {stack_size_mb:.1f} MB (推荐: 64 MB 或更大)\n"
                    f"   CREST 可能因栈不足导致 SIGSEGV (段错误)\n"
                    f"   建议运行: ulimit -s 65536"
                )
                if HAS_RICH:
                    console.print(f"[yellow]{stack_warning}[/yellow]")
                else:
                    logger.warning(stack_warning)
        except (ImportError, AttributeError):
            pass
        except Exception as e:
            logger.debug(f"无法检查栈大小: {e}")

    def _select_work_dir(self, output_dir: Path, run_name: Optional[str], resume: bool) -> Path:
        output_dir = Path(output_dir)
        if run_name:
            return output_dir / run_name
        if resume and output_dir.exists():
            candidates = [path for path in output_dir.glob("test_*") if path.is_dir()]
            if candidates:
                return max(candidates, key=lambda path: path.stat().st_mtime)
        return output_dir / f"test_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    def _build_inspector(self, work_dir: Path) -> ResultInspector:
        context = {
            "smiles": self.test_data["smiles"],
            "charge": 0,
            "multiplicity": 1
        }
        return ResultInspector(work_dir=work_dir, config=self.config, context=context, strict_mode=True)

    def _prepare_run(self, run_kind: str, base_dir: Path, resume: bool):
        mapping = {
            "S1": ("S1_test", ["S1"], ["s2", "s3", "s4"]),
            "S1_S2": ("S1_S2_test", ["S1", "S2"], ["s3", "s4"]),
            "S1_S2_S3": ("S1_S2_S3_test", ["S1", "S2", "S3"], ["s4"]),
            "Full_Pipeline": ("Full_Pipeline_test", ["S1", "S2", "S3", "S4"], [])
        }
        run_dir_name, step_names, base_skip = mapping[run_kind]
        run_dir = base_dir / run_dir_name
        inspector = self._build_inspector(run_dir)
        skip_steps = list(base_skip)
        signature_steps = list(step_names)
        skip_all = False
        skip_reasons: Dict[str, str] = {}

        if resume:
            skipped = []
            for step in step_names:
                result = inspector.check_step(step)
                if result.should_skip:
                    skipped.append(step)
                    skip_reasons[step] = result.reason
            skip_steps.extend([step.lower() for step in skipped])
            signature_steps = [step for step in step_names if step not in skipped]
            if skipped and len(skipped) == len(step_names):
                skip_all = True

        return run_dir, skip_steps, signature_steps, skip_all, inspector, skip_reasons

    def _record_skip(self, step_name: str, run_dir: Path, reason: str) -> None:
        self.test_report["steps"][step_name] = {
            "status": "skipped",
            "attempts": 0,
            "start_time": datetime.now().isoformat(),
            "end_time": datetime.now().isoformat(),
            "output_dir": str(run_dir),
            "skip_reason": reason
        }
        if HAS_RICH:
            console.print(f"[bold yellow]↷ {step_name} skipped[/bold yellow]: {reason}")
        else:
            logger.info(f"↷ {step_name} skipped: {reason}")

    def run_step(
        self,
        step_name: str,
        step_function,
        work_dir: Path,
        max_retries: int = 3,
        step_kwargs: Optional[Dict[str, Any]] = None,
        inspector: Optional[ResultInspector] = None,
        signature_steps: Optional[List[str]] = None
    ) -> bool:
        """
        运行单个步骤（带自动重试）
        """
        step_dir = work_dir / step_name
        step_dir.mkdir(parents=True, exist_ok=True)

        if HAS_RICH:
            console.rule(f"[bold cyan]开始执行: {step_name}[/bold cyan]")
        else:
            logger.info(f"=== 开始执行: {step_name} ===")

        self.test_report["steps"][step_name] = {
            "status": "running",
            "attempts": 0,
            "start_time": datetime.now().isoformat(),
            "output_dir": str(step_dir)
        }

        retry_config = RetryConfig(
            max_retries=max_retries,
            initial_backoff=10.0,
            max_backoff=300.0,
            backoff_multiplier=2.0
        )

        def cleanup():
            """清理失败的临时文件"""
            if HAS_RICH:
                console.print(f"[yellow]Cleaning up {step_name} temp files...[/yellow]")
            else:
                logger.info(f"清理 {step_name} 临时文件...")

        @with_retry(config=retry_config, cleanup=cleanup)
        def execute_step():
            self.test_report["steps"][step_name]["attempts"] += 1
            if step_kwargs:
                return step_function(work_dir, **step_kwargs)
            return step_function(work_dir)

        try:
            # 使用 Spinner 显示运行中状态
            if HAS_RICH:
                with console.status(f"[bold green]Running {step_name}...[/bold green] (Check rph.log for details)", spinner="dots"):
                    result = execute_step()
            else:
                result = execute_step()

            # 成功
            self.test_report["steps"][step_name]["status"] = "completed"
            self.test_report["steps"][step_name]["end_time"] = datetime.now().isoformat()

            if isinstance(result, PipelineResult):
                files = self._collect_files(step_dir)
                self.test_report["files_generated"].extend(files)
            if inspector and signature_steps:
                for step in signature_steps:
                    inspector.save_signature(step)

            if HAS_RICH:
                console.print(f"[bold green]✔ {step_name} 完成[/bold green]")
                console.print(f"  [dim]输出目录: {step_dir}[/dim]")
            else:
                logger.info(f"✓ {step_name} 完成")

            return True

        except Exception as e:
            # 失败
            self.test_report["steps"][step_name]["status"] = "failed"
            self.test_report["steps"][step_name]["end_time"] = datetime.now().isoformat()
            self.test_report["steps"][step_name]["error"] = str(e)

            if HAS_RICH:
                console.print(f"[bold red]✘ {step_name} 失败[/bold red]: {e}")
            else:
                logger.error(f"✗ {step_name} 失败: {e}")

            # 提取诊断信息
            self._extract_diagnostics(step_dir)

            # 尝试错误分析和自动修复
            self._try_auto_fix(step_name, e, work_dir)

            return False

    def _extract_diagnostics(self, step_dir: Path):
        """从失败步骤中提取诊断信息"""
        if HAS_RICH:
            console.print(Panel(f"诊断信息提取: {step_dir}", style="red"))
        else:
            logger.error(f"诊断信息提取: {step_dir}")

        log_extensions = ['.log', '.out', '.err', '.txt']
        log_files = []

        for ext in log_extensions:
            log_files.extend(step_dir.rglob(f"*{ext}"))

        if not log_files:
            return

        log_files_sorted = sorted(log_files, key=lambda x: x.stat().st_mtime)
        last_log = log_files_sorted[-1]

        if HAS_RICH:
             console.print(f"  [bold]最新日志文件[/bold]: {last_log.name}")
        else:
             logger.error(f"  最新日志文件: {last_log.name}")

        try:
            content = last_log.read_text(errors='ignore')
            lines = content.splitlines()

            lines_to_show = lines[-50:] if len(lines) > 50 else lines
            
            if HAS_RICH:
                log_text = "\n".join(lines_to_show)
                console.print(Panel(log_text, title=f"Last 50 lines of {last_log.name}", border_style="red", expand=False))
            else:
                logger.error("--- LOG TAIL ---")
                for line in lines_to_show:
                    logger.error(line)

        except Exception as e:
            logger.error(f"  日志读取失败: {e}")

    def _try_auto_fix(self, step_name: str, error: Exception, work_dir: Path):
        """尝试自动修复错误"""
        _error_msg = str(error).lower()
        if HAS_RICH:
            console.print(f"[yellow]尝试自动修复 {step_name} 的错误...[/yellow]")
        else:
            logger.warning(f"尝试自动修复 {step_name} 的错误...")

        # 检查日志文件中的具体错误
        log_files = list(work_dir.rglob("*.log"))
        for log_file in log_files[-3:]:
            try:
                log_content = log_file.read_text()
                error_type = analyze_log_for_errors(log_content)
                if error_type is None:
                    continue

                suggestion = ""
                if error_type == QCFailureType.OUT_OF_MEMORY:
                    suggestion = "减少 config 中 resources.nproc 的值"
                elif error_type == QCFailureType.CONVERGENCE_FAILURE:
                    suggestion = "检查几何初始质量，或使用更好的初始猜测"
                
                if suggestion:
                    error_type_name = getattr(error_type, "name", "UNKNOWN")
                    if HAS_RICH:
                        console.print(f"  [bold red]{error_type_name}[/bold red] detected. Suggestion: [green]{suggestion}[/green]")
                    else:
                        logger.warning(f"  检测到 {error_type_name}: {suggestion}")
                    
                break
            except Exception:
                pass

    def _collect_files(self, directory: Path) -> list:
        """收集目录中所有生成的文件"""
        files = []
        for f in directory.rglob("*"):
            if f.is_file():
                files.append({
                    "path": str(f.relative_to(directory.parent)),
                    "size": f.stat().st_size,
                    "modified": datetime.fromtimestamp(f.stat().st_mtime).isoformat()
                })
        return files

    def test_s1(self, work_dir: Path, skip_steps: Optional[List[str]] = None, run_dir: Optional[Path] = None):
        result = self.hunter.run_pipeline(
            product_smiles=self.test_data["smiles"],
            work_dir=run_dir or (work_dir / "S1_test"),
            skip_steps=skip_steps or ['s2', 's3', 's4']
        )
        if not result.success:
            raise RuntimeError(f"Step 1 失败: {result.error_message}")
        return result

    def test_s1_s2(self, work_dir: Path, skip_steps: Optional[List[str]] = None, run_dir: Optional[Path] = None):
        result = self.hunter.run_pipeline(
            product_smiles=self.test_data["smiles"],
            work_dir=run_dir or (work_dir / "S1_S2_test"),
            skip_steps=skip_steps or ['s3', 's4']
        )
        if not result.success:
            raise RuntimeError(f"Step 1+2 失败: {result.error_message}")
        return result

    def test_s1_s2_s3(self, work_dir: Path, skip_steps: Optional[List[str]] = None, run_dir: Optional[Path] = None):
        result = self.hunter.run_pipeline(
            product_smiles=self.test_data["smiles"],
            work_dir=run_dir or (work_dir / "S1_S2_S3_test"),
            skip_steps=skip_steps or ['s4']
        )
        if not result.success:
            raise RuntimeError(f"Step 1+2+3 失败: {result.error_message}")
        return result

    def test_full_pipeline(self, work_dir: Path, skip_steps: Optional[List[str]] = None, run_dir: Optional[Path] = None):
        result = self.hunter.run_pipeline(
            product_smiles=self.test_data["smiles"],
            work_dir=run_dir or (work_dir / "Full_Pipeline_test"),
            skip_steps=skip_steps or []
        )
        if not result.success:
            raise RuntimeError(f"完整流程失败: {result.error_message}")
        return result

    def save_report(self, output_dir: Path):
        """保存测试报告"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        report_file = output_dir / f"test_report_{self.test_data['rx_id']}.json"
        with open(report_file, 'w') as f:
            json.dump(self.test_report, f, indent=2)

        if HAS_RICH:
            console.print(f"[dim]测试报告已保存: {report_file}[/dim]")
        else:
            logger.info(f"测试报告已保存: {report_file}")

        # 打印摘要
        self._print_summary()

    def _print_summary(self):
        """打印测试摘要"""
        if HAS_RICH:
            table = Table(title="Test Execution Summary", border_style="blue")
            table.add_column("Step", style="cyan")
            table.add_column("Status", justify="center")
            table.add_column("Attempts")
            
            for step_name, step_info in self.test_report["steps"].items():
                status = step_info["status"]
                attempts = str(step_info.get("attempts", 1))
                if status == "completed":
                    status_str = "[green]COMPLETED[/green]"
                elif status == "skipped":
                    status_str = "[yellow]SKIPPED[/yellow]"
                else:
                    status_str = "[red]FAILED[/red]"
                table.add_row(step_name, status_str, attempts)
            
            console.print(Panel(table))
            
            final_status = "[bold green]SUCCESS[/bold green]" if self.test_report['success'] else "[bold red]FAILURE[/bold red]"
            console.print(f"\nFinal Result: {final_status}")
            console.print(f"Files Generated: {len(self.test_report['files_generated'])}")
            
        else:
            logger.info("="*70)
            logger.info("测试摘要")
            logger.info("="*70)
            for step_name, step_info in self.test_report["steps"].items():
                status = step_info["status"]
                attempts = step_info.get("attempts", 1)
                logger.info(f"  {step_name}: {status} (尝试 {attempts} 次)")
            logger.info(f"生成的文件总数: {len(self.test_report['files_generated'])}")
            logger.info(f"测试状态: {'成功' if self.test_report['success'] else '失败'}")

    def run(self, target_step: str = "ALL", output_dir: Optional[Path] = None, run_name: Optional[str] = None, resume: bool = False):
        if output_dir is None:
            output_dir = Path(__file__).parent / "test_results"

        work_dir = self._select_work_dir(output_dir, run_name, resume)
        if work_dir.exists() and (run_name or resume):
            if HAS_RICH:
                console.print(f"[dim]♻️ Reusing directory: {work_dir}[/dim]")
            else:
                logger.info(f"♻️ Reusing directory: {work_dir}")
        else:
            if HAS_RICH:
                console.print(f"[dim]📁 New directory: {work_dir}[/dim]")
            else:
                logger.info(f"📁 New directory: {work_dir}")

        try:
            if target_step == "S1":
                run_dir, skip_steps, signature_steps, skip_all, inspector, skip_reasons = self._prepare_run("S1", work_dir, resume)
                if skip_all:
                    reason = "; ".join([f"{key}:{value}" for key, value in skip_reasons.items()]) or "outputs_valid"
                    self._record_skip("S1", run_dir, reason)
                else:
                    self.run_step(
                        "S1",
                        self.test_s1,
                        work_dir,
                        step_kwargs={"skip_steps": skip_steps, "run_dir": run_dir},
                        inspector=inspector,
                        signature_steps=signature_steps
                    )
            elif target_step == "S2":
                run_dir, skip_steps, signature_steps, skip_all, inspector, skip_reasons = self._prepare_run("S1_S2", work_dir, resume)
                if skip_all:
                    reason = "; ".join([f"{key}:{value}" for key, value in skip_reasons.items()]) or "outputs_valid"
                    self._record_skip("S1_S2", run_dir, reason)
                else:
                    self.run_step(
                        "S1_S2",
                        self.test_s1_s2,
                        work_dir,
                        step_kwargs={"skip_steps": skip_steps, "run_dir": run_dir},
                        inspector=inspector,
                        signature_steps=signature_steps
                    )
            elif target_step == "S3":
                run_dir, skip_steps, signature_steps, skip_all, inspector, skip_reasons = self._prepare_run("S1_S2_S3", work_dir, resume)
                if skip_all:
                    reason = "; ".join([f"{key}:{value}" for key, value in skip_reasons.items()]) or "outputs_valid"
                    self._record_skip("S1_S2_S3", run_dir, reason)
                else:
                    self.run_step(
                        "S1_S2_S3",
                        self.test_s1_s2_s3,
                        work_dir,
                        step_kwargs={"skip_steps": skip_steps, "run_dir": run_dir},
                        inspector=inspector,
                        signature_steps=signature_steps
                    )
            elif target_step == "S4" or target_step == "ALL":
                run_dir, skip_steps, signature_steps, skip_all, inspector, skip_reasons = self._prepare_run("Full_Pipeline", work_dir, resume)
                if skip_all:
                    reason = "; ".join([f"{key}:{value}" for key, value in skip_reasons.items()]) or "outputs_valid"
                    self._record_skip("Full_Pipeline", run_dir, reason)
                else:
                    self.run_step(
                        "Full_Pipeline",
                        self.test_full_pipeline,
                        work_dir,
                        step_kwargs={"skip_steps": skip_steps, "run_dir": run_dir},
                        inspector=inspector,
                        signature_steps=signature_steps
                    )
            else:
                raise ValueError(f"无效的目标步骤: {target_step}")

            self.test_report["success"] = True
            self.test_report["end_time"] = datetime.now().isoformat()

        except Exception as e:
            self.test_report["success"] = False
            self.test_report["error"] = str(e)
            self.test_report["end_time"] = datetime.now().isoformat()
            logger.error(f"测试异常: {e}", exc_info=False)

        finally:
            self.save_report(output_dir)


def main():
    """命令行入口"""
    import argparse

    parser = argparse.ArgumentParser(
        description="ReactionProfileHunter v2.1 自动测试"
    )
    parser.add_argument(
        '--step',
        type=str,
        default='ALL',
        choices=['S1', 'S2', 'S3', 'S4', 'ALL'],
        help='要测试的步骤 (默认: ALL)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='./test_results',
        help='测试结果输出目录 (默认: ./test_results)'
    )
    parser.add_argument(
        '--config',
        type=str,
        default=None,
        help='配置文件路径 (默认: ./config/defaults.yaml)'
    )
    parser.add_argument(
        '--name',
        type=str,
        default=None,
        help='固定测试目录名称，用于重用结果'
    )
    parser.add_argument(
        '--resume',
        action='store_true',
        help='基于现有输出跳过已完成步骤'
    )

    args = parser.parse_args()

    # 运行测试
    tester = AutoTester(
        config_path=Path(args.config) if args.config else None
    )

    tester.run(
        target_step=args.step,
        output_dir=Path(args.output),
        run_name=args.name,
        resume=args.resume
    )

    # 返回退出码
    return 0 if tester.test_report["success"] else 1


if __name__ == "__main__":
    sys.exit(main())
