"""
parameter_file.py -- reads src/parameters.py as the source of the GUI's forms: every constant is a Parameter
(its default, the inline comment and the indented continuation comments as help text, the section header
'# --- Name [solver|compare|both]: sub-heading ---' as group, mode tag and sub-heading).
"""

import ast
import io
import re
import tokenize
from dataclasses import dataclass
from pathlib import Path

PARAMETERS_FILE = Path(__file__).resolve().parent.parent / "src" / "parameters.py"
SECTION = re.compile(r"^# --- (.*?) -*$")
TAG = re.compile(r"\s*\[(solver|compare|both)\]")


@dataclass(frozen=True)
class Parameter:
    group: str          # expander title ('Run', 'Cooler', ...)
    tag: str            # mode: 'solver', 'compare' or 'both'
    section: str        # sub-heading inside the group
    name: str           # the constant's name in parameters.py
    default: object
    help: str


def read_parameters(path=PARAMETERS_FILE) -> list:
    """The Parameters of parameters.py in file order."""
    source = Path(path).read_text()
    lines = source.splitlines()
    comments = {tok.start[0]: tok.string.lstrip("# ").strip()
                for tok in tokenize.generate_tokens(io.StringIO(source).readline) if tok.type == tokenize.COMMENT}
    entries = []
    for node in ast.parse(source).body:
        if not (isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)):
            continue
        section = next((m.group(1).strip() for l in reversed(lines[:node.lineno - 1]) if (m := SECTION.match(l))), "Other")
        head, _, sub = section.partition(":")
        tag = (TAG.search(head) or [None, "both"])[1]
        help_parts = [comments.get(node.lineno, "")]               # inline comment + indented continuation lines
        n = node.end_lineno + 1
        while n <= len(lines) and lines[n - 1].startswith(" ") and lines[n - 1].strip().startswith("#"):
            help_parts.append(comments[n])
            n += 1
        entries.append(Parameter(group=TAG.sub("", head).strip(), tag=tag, section=sub.strip(), name=node.targets[0].id,
                                 default=ast.literal_eval(node.value), help=" ".join(p for p in help_parts if p)))
    return entries
