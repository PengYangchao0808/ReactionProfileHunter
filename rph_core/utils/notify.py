import logging
import platform
import subprocess
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)


def notify_completion(title: str, message: str, config: Optional[Dict[str, Any]] = None) -> bool:
    notify_cfg = (config or {}).get("run", {}).get("notify", {})
    enabled = notify_cfg.get("enabled", False)
    if not enabled:
        return False
    backend = notify_cfg.get("backend", "auto")
    if backend == "off":
        return False
    if backend in ("auto", "notifypy"):
        if _notify_notifypy(title, message):
            return True
        if backend == "notifypy":
            return False
    if backend in ("auto", "native"):
        if _notify_native(title, message):
            return True
        if backend == "native":
            return False
    logger.info(f"Notification: {title} - {message}")
    return False


def _notify_notifypy(title: str, message: str) -> bool:
    try:
        from notifypy import Notify
    except Exception:
        return False
    try:
        notification = Notify()
        notification.title = title
        notification.message = message
        notification.send(block=False)
        return True
    except Exception:
        return False


def _notify_native(title: str, message: str) -> bool:
    system = platform.system()
    try:
        if system == "Linux":
            subprocess.run(["notify-send", title, message], check=False, timeout=5)
            return True
        if system == "Darwin":
            script = f'display notification "{_escape_applescript(message)}" with title "{_escape_applescript(title)}"'
            subprocess.run(["osascript", "-e", script], check=False, timeout=5)
            return True
        if system == "Windows":
            ps_script = _build_windows_toast(title, message)
            subprocess.run(["powershell", "-NoProfile", "-Command", ps_script], check=False, timeout=8)
            return True
    except Exception:
        return False
    return False


def _escape_applescript(text: str) -> str:
    return text.replace("\\", "\\\\").replace('"', '\\"')


def _build_windows_toast(title: str, message: str) -> str:
    safe_title = title.replace("'", "''")
    safe_message = message.replace("'", "''")
    return (
        "[Windows.UI.Notifications.ToastNotificationManager, Windows.UI.Notifications, ContentType = WindowsRuntime] | Out-Null;"
        "$template = [Windows.UI.Notifications.ToastTemplateType]::ToastText02;"
        "$xml = [Windows.UI.Notifications.ToastNotificationManager]::GetTemplateContent($template);"
        "$text = $xml.GetElementsByTagName('text');"
        f"$text.Item(0).AppendChild($xml.CreateTextNode('{safe_title}')) | Out-Null;"
        f"$text.Item(1).AppendChild($xml.CreateTextNode('{safe_message}')) | Out-Null;"
        "$toast = [Windows.UI.Notifications.ToastNotification]::new($xml);"
        "$notifier = [Windows.UI.Notifications.ToastNotificationManager]::CreateToastNotifier('RPH');"
        "$notifier.Show($toast);"
    )
