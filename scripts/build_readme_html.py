#!/usr/bin/env python3
"""
Markdown -> HTML converter for RPH README files.

Self-contained: stdlib only. Produces standalone .html files in docs/.
Handles: headings, paragraphs, lists (ul/ol), tables (GFM), code blocks,
inline code, bold, italic, links, images, blockquotes, hr, and html-escapes
all content. GitHub-flavored markdown table support is built in.
"""

import re
import sys
import html as html_lib
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DOCS = REPO / "docs"
README = REPO / "README.md"
README_CN = REPO / "README.zh-CN.md"

CSS = """
:root {
  --fg: #1f2328;
  --fg-dim: #57606a;
  --bg: #ffffff;
  --border: #d0d7de;
  --accent: #0969da;
  --code-bg: #f6f8fa;
  --blockquote: #d0d7de;
  --table-stripe: #f6f8fa;
  --badge-bg: #eaeef2;
}
@media (prefers-color-scheme: dark) {
  :root {
    --fg: #e6edf3;
    --fg-dim: #8b949e;
    --bg: #0d1117;
    --border: #30363d;
    --accent: #58a6ff;
    --code-bg: #161b22;
    --blockquote: #30363d;
    --table-stripe: #161b22;
    --badge-bg: #21262d;
  }
}
* { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Noto Sans",
    Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji";
  background: var(--bg);
  color: var(--fg);
  line-height: 1.6;
  margin: 0;
  padding: 2rem 1rem;
  font-size: 16px;
}
.container { max-width: 920px; margin: 0 auto; }
h1, h2, h3, h4, h5, h6 {
  font-weight: 600;
  line-height: 1.25;
  margin-top: 24px;
  margin-bottom: 16px;
}
h1 { font-size: 2em; border-bottom: 1px solid var(--border); padding-bottom: .3em; }
h2 { font-size: 1.5em; border-bottom: 1px solid var(--border); padding-bottom: .3em; }
h3 { font-size: 1.25em; }
h4 { font-size: 1em; }
p { margin: 0 0 16px; }
a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }
code {
  font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas,
    "Liberation Mono", monospace;
  background: var(--code-bg);
  padding: 0.2em 0.4em;
  border-radius: 6px;
  font-size: 85%;
}
pre {
  background: var(--code-bg);
  padding: 16px;
  border-radius: 6px;
  overflow: auto;
  font-size: 85%;
  line-height: 1.45;
}
pre code { background: transparent; padding: 0; font-size: 100%; }
blockquote {
  border-left: 4px solid var(--blockquote);
  color: var(--fg-dim);
  padding: 0 1em;
  margin: 0 0 16px;
}
table {
  border-collapse: collapse;
  width: 100%;
  margin: 0 0 16px;
  display: block;
  overflow: auto;
}
th, td {
  border: 1px solid var(--border);
  padding: 6px 13px;
  text-align: left;
}
tr:nth-child(2n) { background: var(--table-stripe); }
th { font-weight: 600; }
hr { border: 0; border-top: 1px solid var(--border); margin: 24px 0; }
ul, ol { padding-left: 2em; margin: 0 0 16px; }
li + li { margin-top: 0.25em; }
img { max-width: 100%; height: auto; }
.badge {
  display: inline-block;
  padding: 2px 8px;
  border-radius: 12px;
  background: var(--badge-bg);
  font-size: 75%;
  color: var(--fg-dim);
  margin: 0 2px;
}
.center { text-align: center; }
.toc { background: var(--code-bg); padding: 16px; border-radius: 6px; margin: 16px 0; }
.toc ul { margin: 0; }
"""


def esc(s: str) -> str:
    return html_lib.escape(s, quote=True)


def render_inline(text: str) -> str:
    """Render inline markdown: bold, italic, code, links, images, badges."""
    # 1. Extract code spans first (protect from further processing).
    code_spans: list[str] = []

    def _save(m: re.Match) -> str:
        code_spans.append(m.group(1))
        return f"\x00CODE{len(code_spans) - 1}\x00"

    text = re.sub(r"`([^`]+)`", _save, text)

    # 2. Images: ![alt](url)
    text = re.sub(
        r"!\[([^\]]*)\]\(([^)\s]+)(?:\s+\"([^\"]*)\")?\)",
        lambda m: f'<img src="{esc(m.group(2))}" alt="{esc(m.group(1))}"' +
        (f' title="{esc(m.group(3))}"' if m.group(3) else "") + ">",
        text,
    )

    # 3. Links: [text](url)
    text = re.sub(
        r"\[([^\]]+)\]\(([^)\s]+)(?:\s+\"([^\"]*)\")?\)",
        lambda m: f'<a href="{esc(m.group(2))}"' +
        (f' title="{esc(m.group(3))}"' if m.group(3) else "") +
        f">{m.group(1)}</a>",
        text,
    )

    # 4. Bold (**...** or __...__)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"__(.+?)__", r"<strong>\1</strong>", text)

    # 5. Italic (*...* or _..._) - careful not to break bold markers
    text = re.sub(r"(?<!\*)\*([^\*\n]+?)\*(?!\*)", r"<em>\1</em>", text)
    text = re.sub(r"(?<!_)_([^_\n]+?)_(?!_)", r"<em>\1</em>", text)

    # 6. Strikethrough ~~...~~
    text = re.sub(r"~~(.+?)~~", r"<del>\1</del>", text)

    # 7. Restore code spans
    for i, c in enumerate(code_spans):
        text = text.replace(f"\x00CODE{i}\x00", f"<code>{esc(c)}</code>")

    return text


def render_table(lines: list[str]) -> str:
    """Render a GFM table block. Caller has already stripped the block."""
    if len(lines) < 2:
        return ""
    header_cells = [c.strip() for c in lines[0].strip().strip("|").split("|")]
    # Skip alignment row (lines[1])
    out = ["<table>"]
    out.append("<thead><tr>")
    for h in header_cells:
        out.append(f"<th>{render_inline(h)}</th>")
    out.append("</tr></thead><tbody>")
    for line in lines[2:]:
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        out.append("<tr>")
        for c in cells:
            out.append(f"<td>{render_inline(c)}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "\n".join(out)


def is_table_separator(line: str) -> bool:
    s = line.strip().strip("|")
    if not s:
        return False
    cells = [c.strip() for c in s.split("|")]
    return all(re.fullmatch(r":?-{3,}:?", c) for c in cells)


def md_to_html(md: str) -> str:
    """Convert markdown source to HTML body."""
    lines = md.splitlines()
    out: list[str] = []
    i = 0
    in_html_block = False

    def flush_para(buf: list[str]) -> None:
        if not buf:
            return
        text = " ".join(buf).strip()
        if text:
            out.append(f"<p>{render_inline(text)}</p>")
        buf.clear()

    para_buf: list[str] = []
    list_stack: list[str] = []  # 'ul' or 'ol'

    def close_lists() -> None:
        while list_stack:
            out.append(f"</{list_stack.pop()}>")

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Fenced code block
        if stripped.startswith("```"):
            flush_para(para_buf)
            close_lists()
            lang = stripped[3:].strip()
            i += 1
            code_lines: list[str] = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code_lines.append(lines[i])
                i += 1
            i += 1  # skip closing ```
            code_text = "\n".join(code_lines)
            klass = f' class="language-{esc(lang)}"' if lang else ""
            out.append(f"<pre><code{klass}>{esc(code_text)}</code></pre>")
            continue

        # HTML block (raw passthrough until blank)
        if stripped.startswith("<") and not stripped.startswith("<!--"):
            flush_para(para_buf)
            close_lists()
            html_lines = [line]
            i += 1
            while i < len(lines) and lines[i].strip() and (
                lines[i].lstrip().startswith("<") or lines[i].strip().startswith("<")
            ):
                html_lines.append(lines[i])
                i += 1
            out.append("\n".join(html_lines))
            continue

        # Heading
        m = re.match(r"^(#{1,6})\s+(.*)", line)
        if m:
            flush_para(para_buf)
            close_lists()
            level = len(m.group(1))
            content = m.group(2).strip()
            # Auto-generate id from text
            anchor = re.sub(r"[^\w\s-]", "", content.lower())
            anchor = re.sub(r"\s+", "-", anchor)
            out.append(f'<h{level} id="{esc(anchor)}">{render_inline(content)}</h{level}>')
            i += 1
            continue

        # Horizontal rule
        if re.fullmatch(r"\s*([-*_])\s*\1\s*\1[\s\1]*", line):
            flush_para(para_buf)
            close_lists()
            out.append("<hr>")
            i += 1
            continue

        # Blockquote
        if stripped.startswith(">"):
            flush_para(para_buf)
            close_lists()
            quote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                quote_lines.append(re.sub(r"^>\s?", "", lines[i].lstrip()))
                i += 1
            inner = md_to_html("\n".join(quote_lines))
            out.append(f"<blockquote>\n{inner}\n</blockquote>")
            continue

        # GFM Table
        if "|" in line and i + 1 < len(lines) and is_table_separator(lines[i + 1]):
            flush_para(para_buf)
            close_lists()
            tbl = [line]
            i += 1
            tbl.append(lines[i])  # separator
            i += 1
            while i < len(lines) and "|" in lines[i] and lines[i].strip():
                tbl.append(lines[i])
                i += 1
            out.append(render_table(tbl))
            continue

        # Lists
        ol_m = re.match(r"^(\s*)(\d+)\.\s+(.*)", line)
        ul_m = re.match(r"^(\s*)[-*+]\s+(.*)", line)
        if ol_m or ul_m:
            flush_para(para_buf)
            indent, content = (ol_m.group(1), ol_m.group(3)) if ol_m else (
                ul_m.group(1), ul_m.group(2)
            )
            kind = "ol" if ol_m else "ul"
            depth = len(indent) // 2
            # Open/close lists to match depth and kind
            while len(list_stack) > depth + 1:
                out.append(f"</{list_stack.pop()}>")
            if len(list_stack) <= depth:
                while len(list_stack) < depth:
                    out.append("<ul>")
                    list_stack.append("ul")
                if not list_stack or list_stack[-1] != kind:
                    if list_stack:
                        out.append(f"</{list_stack.pop()}>")
                    out.append(f"<{kind}>")
                    list_stack.append(kind)
            out.append(f"<li>{render_inline(content)}</li>")
            i += 1
            continue

        # Blank line
        if not stripped:
            flush_para(para_buf)
            close_lists()
            i += 1
            continue

        # Paragraph (accumulate)
        para_buf.append(line)
        i += 1

    flush_para(para_buf)
    close_lists()
    return "\n".join(out)


def build_page(title: str, body_html: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{esc(title)}</title>
<style>{CSS}</style>
</head>
<body>
<div class="container">
{body_html}
</div>
</body>
</html>
"""


def convert(src: Path, dst: Path, title: str) -> None:
    md = src.read_text(encoding="utf-8")
    body = md_to_html(md)
    page = build_page(title, body)
    dst.write_text(page, encoding="utf-8")
    print(f"  {src.name} -> {dst.relative_to(REPO)} ({len(page):,} bytes)")


def main() -> int:
    if not DOCS.exists():
        DOCS.mkdir(parents=True, exist_ok=True)

    print("Converting README to HTML:")
    convert(README, DOCS / "README.html", "ReactionProfileHunter — Documentation")
    convert(README_CN, DOCS / "README.zh-CN.html", "ReactionProfileHunter — 文档")
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
