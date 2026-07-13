# RPH V4 CLI UI 美化方案（修订版）

> 目标：将 V4 的命令行界面从当前平铺的 Python logging 输出，恢复到 V3 风格的 Rich 终端体验（彩色阶段头、进度表格、图标语义、去重刷新），同时保留 `rph_v4.log` 的纯文本格式与 `events.jsonl`/`status.json` 的持久化数据。

---

## 1. 当前问题诊断

### 1.1 V4 实际输出样式

当前执行 `bin/rph_run --csv ... --rx-id 1 --output /tmp/rph4_rx_1 --stop-after s3` 时，终端输出为：

```text
2026-07-12 21:37:42,691 | INFO | rph_core.v4_orchestrator | [V4] Pipeline log: /tmp/rph4_rx_1/rph_v4.log
2026-07-12 21:37:42,716 | INFO | rph_core.utils.config_loader | 已加载配置文件: .../config/defaults.yaml
2026-07-12 21:37:42,729 | INFO | rph_core.utils.resource_utils | Found: crest (config): /opt/software/crest/crest
2026-07-12 21:37:42,756 | INFO | rph_core.utils.resource_utils | Found: xtb (config): /opt/software/xtb/bin/xtb
```

### 1.2 已发现的缺陷

- 当前 `rph_core/v4_orchestrator.py` 的 `_configure_v4_logging()` 使用 `logging.basicConfig()` 输出平铺文本，没有使用 Rich 控制台 handler。`rph_v4.log` 的 `FileHandler` 已正确注册到 `root` logger（[line 1062](rph_core/v4_orchestrator.py:1062)），所以**日志文件落盘是正常的**；之前版本对未注册的判断已过时，不应再“修复”该 bug，避免重复注册或破坏现有 handler。
- `rph_core/utils/log_manager.py` 已封装 `setup_logger()`，但 V4 入口未使用；`setup_logger()` 的实现若直接 `logger.handlers.clear()` 会清掉 pytest/宿主预置的 handler，需要改造为只管理自身创建的 handler。
- V4 各阶段的数据层已就绪，但日志/控制台 UI 没有消费这些结构化事件。

### 1.3 可复用的 V3 资产（仅提取通用资产，不复用业务逻辑）

- `rph_core/utils/v3_progress.py`：可借鉴 `V3StageState` / `V3ItemState` 的“图标 + 颜色 + 标签”设计模式，但**不直接复用这些类**，避免引入 V3 编排语义。
- `rph_core/utils/v3_stage_display.py`：可参考其 Rich 渲染、去重缓存、表格布局的思路，但**不继承 `V3StageDisplay` 或其业务方法**。
- `rph_core/utils/shared_console.py`：Rich `Console` 单例和 `RPH_THEME` 可直接复用。
- `requirements.txt` 已包含 `rich>=13.0`。

### 1.4 V4 数据层已就绪

- `rph_core/utils/stage_progress.py`：S0–S3 的 `StageProgressReporter`，写 `events.jsonl` + `status.json`。`status.json` 的 `structures` 是 **dict**，键为 `structure_id`。
- `rph_core/utils/s4_progress.py`：S4 的 `S4ProgressReporter`，`status.json` 的 `structures` 是 **list**，元素为 dict。
- `rph_core/steps/step4_highlevel/engine.py`：S4 内部创建 `S4ProgressReporter`，编排器只调用 `HighLevelEngine.run()`（[v4_orchestrator.py:648](rph_core/v4_orchestrator.py:648)），没有直接的事件通道。
- `rph_core/v4_watch.py`：读取状态文件做实时展示，但使用原始 ANSI 字符串；`--overview` 分支在 `--watch` 时不会检查终态，可能无限循环（[v4_watch.py:226](rph_core/v4_watch.py:226)）。

### 1.5 状态模型不统一

当前各阶段/任务出现的状态字符串至少包括：

| 状态 | 出现位置 |
|------|----------|
| `pending` | 通用 |
| `running` | 通用 |
| `complete` | S3/S4 结构、任务 |
| `completed` | S4 stage 状态 |
| `completed_with_failures` | S3/S4 stage 状态 |
| `failed` | 通用 |
| `degraded` | S4 结构 |
| `opt_failed_sp_complete` | S4 结构 |
| `ts_frequency_unverified` | S4 结构 |
| `skipped` | 任务 |
| `cached` / `reused` | 复用场景 |

UI 层必须增加**统一状态适配器**，将这些字符串映射到稳定的显示状态（pending / running / complete / failed / degraded / cached），再映射到图标和颜色。

---

## 2. 设计目标

| 目标 | 说明 |
|------|------|
| 彩色阶段头 | S0–S4 每个阶段启动/完成有醒目的 Panel/标题 |
| 图标语义 | `→` 运行、`✓` 完成、`✗` 失败、`↻` 复用、`⚠` 降级 |
| 结构表格 | S3/S4 显示 id、kind、source、opt/freq/sp 状态、能量 |
| 不刷屏 | 结构/任务更新时只做增量/去重输出 |
| 日志兼容 | `rph_v4.log` 保持纯文本，方便 grep；log 中无 ANSI/markup |
| 环境兼容 | 非 TTY / CI 自动降级为 plain 输出；无 Rich 时也能运行 |
| Live 监控 | `rph_watch` 升级为 Rich 实时 dashboard，并能在所有阶段完成后自动退出 |
| 统一状态 | S0–S4 的结构/任务状态经 adapter 后使用同一模型渲染 |

---

## 3. 设计系统

### 3.1 配色扩展

复用并扩展 `rph_core/utils/shared_console.py` 的 `RPH_THEME`：

```python
RPH_THEME = Theme({
    "info": "dim cyan",
    "warning": "yellow",
    "error": "bold red",
    "success": "bold green",
    "step.header": "bold magenta",
    "step.title": "bold white",
    "status.pending": "dim",
    "status.running": "cyan",
    "status.complete": "green",
    "status.failed": "red",
    "status.cached": "blue",
    "status.degraded": "yellow",
    "task.opt": "bright_blue",
    "task.freq": "bright_magenta",
    "task.sp": "bright_cyan",
    "energy": "green",
})
```

### 3.2 统一状态与图标

新建 `rph_core/utils/ui_state.py`：

```python
from enum import Enum

class UiStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"
    DEGRADED = "degraded"
    CACHED = "cached"

_STATUS_MAP: dict[str, UiStatus] = {
    "pending": UiStatus.PENDING,
    "running": UiStatus.RUNNING,
    "complete": UiStatus.COMPLETE,
    "completed": UiStatus.COMPLETE,
    "completed_with_failures": UiStatus.COMPLETE,
    "failed": UiStatus.FAILED,
    "degraded": UiStatus.DEGRADED,
    "opt_failed_sp_complete": UiStatus.DEGRADED,
    "ts_frequency_unverified": UiStatus.DEGRADED,
    "cached": UiStatus.CACHED,
    "reused": UiStatus.CACHED,
    "skipped": UiStatus.CACHED,
}

_STATUS_STYLE: dict[UiStatus, tuple[str, str]] = {
    UiStatus.PENDING:  ("○", "dim"),
    UiStatus.RUNNING:  ("→", "cyan"),
    UiStatus.COMPLETE: ("✓", "green"),
    UiStatus.FAILED:   ("✗", "red"),
    UiStatus.DEGRADED: ("⚠", "yellow"),
    UiStatus.CACHED:   ("↻", "blue"),
}

def normalize_status(raw: str | None) -> UiStatus:
    return _STATUS_MAP.get((raw or "").strip().lower(), UiStatus.PENDING)

def status_markup(status: UiStatus) -> str:
    icon, color = _STATUS_STYLE[status]
    return f"[{color}]{icon}[/]"
```

---

## 4. 日志与控制台输出改造

### 4.1 设计原则

- **logger 层级**：配置 `root` logger（或 `rph_core` logger），因为 V4 代码使用 `logging.getLogger(__name__)`（如 `rph_core.v4_orchestrator`），都在 `rph_core` 层级下。
- **不破坏宿主 handler**：不清空已有 handlers；只添加/移除自己创建的 handler，并通过 `handler._rph_v4_handler = True` 标记去重。
- **文件 handler 去重**：保留现有 `_rph_v4_log_path` 检查逻辑，避免同一进程重复注册同一路径的 `FileHandler`。
- **Rich 降级**：非 TTY 或无 `rich` 时，控制台使用 plain `StreamHandler`；文件 handler 始终为 plain。

### 4.2 入口改动

`rph_core/v4_orchestrator.py` 的 `main()` 调用新的 `setup_v4_logging()`：

```python
from rph_core.utils.log_manager import setup_v4_logging

def main(argv=None):
    parser = argparse.ArgumentParser(...)
    parser.add_argument("--no-color", action="store_true", help="Disable colored console output")
    ...
    args = parser.parse_args(argv)

    log_file = args.log_file or args.output / "rph_v4.log"
    color = not args.no_color and sys.stdout.isatty()
    setup_v4_logging(
        log_file=log_file,
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        rich_console=color,
    )
```

### 4.3 `setup_v4_logging()` 实现

```python
def setup_v4_logging(
    log_file: Path | None,
    level: int = logging.INFO,
    rich_console: bool = True,
) -> None:
    """Configure root logger for V4: plain file + Rich/plain console."""
    root = logging.getLogger()
    root.setLevel(min(level, root.level or logging.WARNING))

    # 1. Console handler: 最多一个 V4 console handler，Rich 或 plain
    for handler in list(root.handlers):
        if getattr(handler, "_rph_v4_console", False):
            root.removeHandler(handler)
            handler.close()

    if rich_console and HAS_RICH:
        console_handler = RichHandler(
            console=get_console(),
            rich_tracebacks=True,
            show_time=True,
            show_path=False,
            markup=True,
        )
    else:
        formatter = logging.Formatter(
            "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
        )
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
    console_handler.setLevel(level)
    console_handler._rph_v4_console = True  # type: ignore[attr-defined]
    root.addHandler(console_handler)

    # 2. File handler: 去重，不删除其他 handler
    if log_file:
        resolved = Path(log_file).resolve()
        for handler in root.handlers:
            if getattr(handler, "_rph_v4_log_path", None) == str(resolved):
                return  # already registered
        resolved.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(resolved, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
        ))
        file_handler._rph_v4_log_path = str(resolved)  # type: ignore[attr-defined]
        root.addHandler(file_handler)
```

### 4.4 废弃 `_configure_v4_logging()`

将 `_configure_v4_logging()` 替换为对 `setup_v4_logging()` 的调用，避免两套 handler 管理逻辑并存。旧函数的去重和文件落盘逻辑已由新函数继承。

---

## 5. 统一数据模型与状态适配

### 5.1 通用 UI 结构模型

新建 `rph_core/utils/ui_adapter.py`：

```python
from dataclasses import dataclass, field
from typing import Any

@dataclass
class UiTask:
    name: str
    status: UiStatus
    engine: str = ""
    method: str = ""
    basis: str = ""
    solvent: str = ""
    output: str = ""
    energy_hartree: float | None = None
    error: str = ""

@dataclass
class UiStructure:
    id: str
    kind: str
    source: str
    status: UiStatus
    current_task: str | None
    tasks: dict[str, UiTask] = field(default_factory=dict)
    usable_for_ml: bool = False
    error: str | None = None
    energy_hartree: float | None = None
```

### 5.2 适配 S3 与 S4 的差异

```python
def adapt_s3_structures(raw: dict[str, Any]) -> list[UiStructure]:
    structures = raw.get("structures", {})
    return [adapt_one(structure_id, data) for structure_id, data in structures.items()]

def adapt_s4_structures(raw: dict[str, Any]) -> list[UiStructure]:
    return [adapt_one(row["id"], row) for row in raw.get("structures", [])]

def adapt_one(structure_id: str, data: dict[str, Any]) -> UiStructure:
    tasks = {}
    for name, task in (data.get("tasks") or {}).items():
        tasks[name] = UiTask(
            name=name,
            status=normalize_status(task.get("status")),
            engine=task.get("engine", ""),
            method=task.get("method", ""),
            basis=task.get("basis", ""),
            energy_hartree=task.get("energy_hartree"),
        )
    return UiStructure(
        id=structure_id,
        kind=data.get("kind", "minimum"),
        source=data.get("input_source", "unknown"),
        status=normalize_status(data.get("status")),
        current_task=data.get("current_task"),
        tasks=tasks,
        usable_for_ml=bool(data.get("usable_for_ml")),
        error=data.get("error"),
        energy_hartree=data.get("sp_energy_hartree") or data.get("energy_hartree"),
    )
```

### 5.3 阶段状态统一

```python
def adapt_stage_status(raw: dict[str, Any]) -> UiStatus:
    return normalize_status(raw.get("status"))
```

---

## 6. RichReporter 与 S0–S4 事件接入

### 6.1 RichReporter 设计

新建 `rph_core/utils/ui_reporter.py`：

```python
class RichReporter:
    """Single data entry for V4 terminal UI.

    Accepts normalized events from StageProgressReporter (S0–S3) and
    S4ProgressReporter (S4) via callbacks, and renders Rich panels/tables.
    """

    def __init__(self, console: Console, color: bool = True, quiet: bool = False):
        self.console = console
        self.color = color
        self.quiet = quiet
        self._stage_cache: dict[str, str] = {}
        self._structure_cache: dict[str, str] = {}

    def stage_started(self, stage: str, meta: dict[str, Any]) -> None:
        ...

    def stage_finished(self, stage: str, status: UiStatus, summary: dict[str, Any]) -> None:
        ...

    def structure_started(self, stage: str, structure: UiStructure) -> None:
        ...

    def structure_updated(self, stage: str, structure: UiStructure) -> None:
        ...

    def structure_finished(self, stage: str, structure: UiStructure) -> None:
        ...
```

### 6.2 S0–S3 接入

`StageProgressReporter` 增加可选的 `event_callback`：

```python
class StageProgressReporter:
    def __init__(..., event_callback: Callable[[str, dict], None] | None = None):
        self._event_callback = event_callback

    def emit(self, event: str, structure_id: str | None = None, **fields):
        ...
        if self._event_callback:
            try:
                self._event_callback(event, record)
            except Exception:
                pass  # UI failures must not abort QC
```

编排器实例化 `StageProgressReporter` 时传入 `RichReporter` 的回调方法。

### 6.3 S4 接入

`HighLevelEngine.run()` 增加可选的 `event_callback`：

```python
class HighLevelEngine:
    def run(
        self,
        structures: Iterable[dict[str, Any]],
        output_dir: Path,
        event_callback: Callable[[str, dict], None] | None = None,
    ) -> Path:
        reporter = S4ProgressReporter(output_dir, materialized, event_callback=event_callback)
        ...
```

`S4ProgressReporter` 在 `emit()` 中同步调用 `event_callback`，与 S0–S3 保持一致。

编排器调用：

```python
s4_manifest = HighLevelEngine(self.config).run(
    s4_structures,
    s4_dir,
    event_callback=rich_reporter.s4_event_callback,
)
```

### 6.4 唯一数据入口原则

- **运行时**：RichReporter 优先通过 `event_callback` 实时消费事件。
- **Fallback / watch**：当无法注入回调（例如外部进程 `rph_watch`）时，读取 `status.json` 并通过 `ui_adapter` 渲染。
- 所有 S0–S4 结构最终都转换为 `UiStructure` 后进入同一渲染路径。

---

## 7. 阶段与结构展示

### 7.1 阶段启动/完成

每个阶段输出一个 Rich `Panel`：

```text
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃  S0 机制验证      RXN_1      [→] running           ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
```

阶段完成：

```text
✓ S1 构象搜索完成    2/2 variants | 0 failed | 1 reused
```

### 7.2 阶段内表格

```text
 variant    smiles                 status     detail
 ─────────────────────────────────────────────────────
 precursor  C1=CC=CC=C1            ✓ complete selected=conf_000
 product    C1CCCCC1               → running  CREST sampling
```

### 7.3 S3/S4 结构表格

```text
 id                 kind    source     status    active task    opt       freq      sp        energy
 ─────────────────────────────────────────────────────────────────────────────────────────────────────
 product            minimum selected   ✓ complete              ✓ done    -         ✓ done   -312.456
 product_int        minimum peb       → running   opt_ts        → running -         -        -
 product_ts         ts      peb        ✓ complete  freq_verify  ✓ done   ✓ done    ✓ done   -312.112
```

---

## 8. 防刷屏机制

参考 `v3_stage_display.py` 的 `_emit_if_changed()`，但作用在 `UiStructure` 指纹上：

```python
def _fingerprint(structure: UiStructure) -> str:
    task_parts = "|".join(
        f"{name}:{task.status.value}:{task.energy_hartree}"
        for name, task in sorted(structure.tasks.items())
    )
    return (
        f"{structure.id}:{structure.status.value}:"
        f"{structure.current_task}:{task_parts}:"
        f"{structure.error}:{structure.energy_hartree}"
    )
```

- `RichReporter` 为每个 `(stage, structure_id)` 缓存最近指纹，只有变化才重新渲染该行。
- 阶段级头只在阶段切换时打印一次。
- 大量 structures 时，默认只显示当前 `running` 的 5 行 + 汇总；`--verbose` 才全量刷新。

---

## 9. `rph_watch` 实时面板

### 9.1 升级为 Rich Live Dashboard

将 `v4_watch.py` 从 ANSI 字符串升级为 `rich.live.Live` + `rich.layout.Layout`：

```text
┌─ Cross-Stage Overview ──┬─ Recent Events ─┐
│ S0 ✓ complete           │ 21:37:42 s1_s.. │
│ S1 ✓ complete           │ 21:37:45 s1_f.. │
│ S2 → running              │                 │
└──────────┬────────────────┴─────────────────┘
┌─ S3 Structures ───────────────────────────────┐
│ id       kind    status    opt    freq  sp  │
└─────────────────────────────────────────────┘
```

### 9.2 修复 `--watch --overview` 终止行为

增加终态检测：

```python
_TERMINAL = {"complete", "completed", "completed_with_failures", "failed"}

def _is_terminal(stage_data: dict[str, Any]) -> bool:
    return normalize_status(stage_data.get("status")).value in _TERMINAL

def _all_terminal(stages: dict[str, Any]) -> bool:
    return all(_is_terminal(data) for data in stages.values())
```

在 `while` 循环中：

```python
while True:
    if args.overview:
        stages = scan_all_stages(args.output)
        print(render_overview(stages, colour), flush=True)
        if args.watch and _all_terminal(stages):
            break
    else:
        status = load_status(args.output)
        ...
    if not args.watch:
        return 0
    time.sleep(args.interval)
```

保留 `--no-color`、`--overview`、`--watch`、`--interval`。

---

## 10. 配置开关

CLI 新增：

| 标志 | 作用 |
|------|------|
| `--no-color` | 强制 plain text |
| `--quiet` / `-q` | 只显示阶段级事件，不显示结构表格 |
| `--verbose` / `-v` | 显示全部结构更新 |
| `--progress {rich,plain,off}` | 进度渲染模式，默认 `rich` |

可选在 `config/defaults.yaml` 中增加：

```yaml
ui:
  theme: default
  color: true
  progress: rich
  quiet: false
```

---

## 11. 实施阶段

### Phase 1：日志基础设施改造（1–2 天）

- 将 `_configure_v4_logging()` 替换为 `setup_v4_logging()`，配置 `root` logger，不破坏已有 handler，支持 Rich/plain 自动切换。
- 新增 `--no-color` 参数。
- 验证 `rph_v4.log` 为纯文本，控制台有颜色，重复初始化不重复注册 handler。

验证：

```bash
python -m py_compile rph_core/v4_orchestrator.py rph_core/utils/log_manager.py
bin/rph_run --csv data/reaxys_cleaned.csv --rx-id 1 --output /tmp/rph4_rx_1 --stop-after s1
# 检查 /tmp/rph4_rx_1/rph_v4.log 是纯文本且无 ANSI 转义序列
grep -P '[\x1b\x9b]' /tmp/rph4_rx_1/rph_v4.log && echo "ANSI found" || echo "clean"
```

### Phase 2：统一状态模型与 UI 渲染层（2–3 天）

- 新建 `rph_core/utils/ui_state.py`：统一状态枚举、图标、颜色。
- 新建 `rph_core/utils/ui_adapter.py`：把 S3 dict/S4 list 差异抹平为 `UiStructure`。
- 新建 `rph_core/utils/ui_reporter.py`：带去重的 `RichReporter`。
- 在 `v4_orchestrator.py` 中实例化 `RichReporter` 并传入 S0–S3 reporter；S4 通过 `event_callback` 接入。

验证：

```bash
bin/rph_run --csv data/reaxys_cleaned.csv --rx-id 1 --output /tmp/rph4_rx_1 --stop-after s3
# 检查终端阶段头、结构表格、去重
```

### Phase 3：改造 `rph_watch`（2 天）

- 用 `rich.live.Live` + `Layout` 重写 `render_status()` / `render_overview()`。
- 使用 `ui_adapter` 统一读取 S0–S4 `status.json`。
- 修复 `--watch --overview` 终态退出。
- 支持 `--no-color` 降级。

验证：

```bash
# 终端 1：完整运行到 S4
bin/rph_run --csv data/reaxys_cleaned.csv --rx-id 1 --output /tmp/rph4_rx_1 --stop-after s4
# 终端 2：watch overview 应自动退出
bin/rph_watch --output /tmp/rph4_rx_1 --watch --overview
```

### Phase 4：测试与文档（1–2 天）

新增 `tests/test_v4_ui.py` 覆盖：

- 非 TTY 下自动降级为 plain
- 无 `rich` 环境仍可用 plain fallback
- `rph_v4.log` 中无 ANSI 转义序列
- 重复调用 `setup_v4_logging()` 不重复注册 handler
- S4 实时事件通过 `event_callback` 到达 `RichReporter`
- `--stop-after s1` 验证（CLI 当前不支持 s0）
- `rph_watch --watch --overview` 在全部阶段终态后自动退出
- `ui_adapter` 对 S3 dict 和 S4 list 输出一致字段

文档：

- 更新 `docs/RPH_V4_COMPLETION_MASTER_PLAN.md` 和 `AGENTS.md` 的 UI 约定。
- 确保 `scripts/ci/check_imports.py rph_core` 通过。

---

## 12. 文件改动清单

| 文件 | 动作 | 说明 |
|------|------|------|
| `rph_core/v4_orchestrator.py` | 修改 | 改用 `setup_v4_logging()`，注入 `RichReporter`，S4 传入 `event_callback` |
| `rph_core/utils/log_manager.py` | 修改 | 新增 `setup_v4_logging()`：root logger、不破坏 host handler、Rich/plain 切换 |
| `rph_core/utils/shared_console.py` | 修改 | 扩展 `RPH_THEME` |
| `rph_core/utils/ui_state.py` | 新建 | 统一状态枚举、图标、颜色、映射 |
| `rph_core/utils/ui_adapter.py` | 新建 | S3/S4 `status.json` → `UiStructure` 适配器 |
| `rph_core/utils/ui_reporter.py` | 新建 | 去重 Rich 渲染器，唯一数据入口 |
| `rph_core/utils/stage_progress.py` | 修改 | 增加 `event_callback` 参数，把事件推给 `RichReporter` |
| `rph_core/utils/s4_progress.py` | 修改 | 增加 `event_callback` 参数，把事件同步给 `RichReporter` |
| `rph_core/steps/step4_highlevel/engine.py` | 修改 | `run()` 接受 `event_callback` 并传给 `S4ProgressReporter` |
| `rph_core/v4_watch.py` | 修改 | 用 Rich live dashboard 重写，修复 `--watch --overview` 终态退出 |
| `config/defaults.yaml` | 可选修改 | 增加 `ui:` 配置块 |
| `tests/test_v4_ui.py` | 新建 | UI 行为测试 |
| `docs/` 与 `AGENTS.md` | 修改 | 更新 UI 约定 |

---

## 13. 风险与注意事项

- **日志文件格式**：`rph_v4.log` 必须保持纯文本，不能带 Rich markup 或 ANSI 转义。`FileHandler` 使用标准 `Formatter`。
- **CI / 非 TTY**：默认应自动降级为 plain；`--no-color` 强制 plain；无 `rich` 时仍可用 plain fallback。
- **Handler 不破坏宿主**：`setup_v4_logging()` 只管理自身标记的 handler，不清空 pytest / Jupyter / 其他框架预置的 handler。
- **去重缓存**：在 `RichReporter` 内维护，不要污染业务对象；`v4_orchestrator` 中避免直接 `print`。
- **性能**：大量结构频繁更新时，关闭 Live 刷新，使用单行 `console.print()` 增量输出，避免重绘整个终端。
- **向后兼容**：`events.jsonl` / `status.json` 的 schema 不变；`event_callback` 为可选参数，不影响现有调用方。
- **状态映射统一**：`ui_state.normalize_status()` 必须覆盖所有已有状态字符串；新增状态字符串时同步更新映射。
- **导入检查**：新增 `rph_core.utils.ui_*` 模块后，确保 `scripts/ci/check_imports.py rph_core` 通过。
- **测试环境**：当前环境中 `rdkit` 缺失会导致 pytest 收集失败；UI 测试应独立或在 `rdkit` 可用环境中运行。
- **不复用 V3 业务类**：只借鉴 V3 的图标/颜色/去重思想；不继承 `V3StageDisplay`、`V3RunProgress` 等包含 V3 编排语义的类。

---

## 14. 预期终端效果（示例）

```text
(rph4) ➜  RPH_V4.0.0 bin/rph_run --csv data/reaxys_cleaned.csv --rx-id 1 --output /tmp/rph4_rx_1 --stop-after s3

[22:15:03] ╔══════════════════════════════════════════════════════╗
           ║  RPH V4 Pipeline   RXN_1   [→] running               ║
           ╚══════════════════════════════════════════════════════╝

[22:15:03] ✓ Loaded config: config/defaults.yaml
[22:15:03] ✓ Found: crest (config) /opt/software/crest/crest
[22:15:03] ✓ Found: xtb (config) /opt/software/xtb/bin/xtb

[22:15:03] ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
           ┃  S0 机制验证      RXN_1      [✓] completed          ┃
           ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
           S0 restored from trusted reaction record rx_id=1

[22:15:04] ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
           ┃  S1 构象搜索      [→] running                          ┃
           ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
           variant    smiles                 status     detail
           ─────────────────────────────────────────────────────
           precursor  C1=CC=CC=C1            ✓ complete selected=conf_000
           product    C1CCCCC1               → running  CREST sampling

[22:18:12] ✓ S1 构象搜索完成    2/2 variants | 0 failed | 1 reused

[22:18:12] ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
           ┃  S2 PEB 扫描      [→] running                          ┃
           ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
           ...

[22:21:45] ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
           ┃  S3 低精度 QC     [→] running                          ┃
           ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
           id                 kind    source     status    active task    opt       freq      sp        energy
           ─────────────────────────────────────────────────────────────────────────────────────────────────────
           product            minimum selected   ✓ complete              ✓ done    -         ✓ done   -312.456
           product_int        minimum peb       → running   opt_ts        → running -         -        -
           product_ts         ts      peb        ✓ complete  freq_verify  ✓ done   ✓ done    ✓ done   -312.112

[22:35:10] ✓ Pipeline stopped after S3
           artifacts: /tmp/rph4_rx_1
           log: /tmp/rph4_rx_1/rph_v4.log
```
