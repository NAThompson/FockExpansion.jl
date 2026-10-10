#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "sympy==1.14.0",
#   "python-flint==0.9.0",
#   "mpmath==1.3.0",
#   "rich==14.1.0",
# ]
# ///
"""Exact recurrence certificate with readable terminal progress and diagnostics."""
from pathlib import Path
import subprocess
import sys
import time
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

HERE = Path(__file__).resolve().parent


def run_stage(script, label, console, log):
    started = time.perf_counter()
    console.print(f"[bold]{label}[/bold]")
    with console.status(label, spinner="dots") as status:
        process = subprocess.Popen(
            [sys.executable, "-u", str(HERE / script)], cwd=HERE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )
        output = []
        for line in process.stdout:
            output.append(line)
            log.write(line)
            log.flush()
            if line.startswith("Running "):
                status.update(line.strip())
                if not console.is_terminal:
                    console.print(Text("  " + line.strip(), style="dim"))
        code = process.wait()
    if code:
        console.print(Panel(Text("".join(output[-40:])), title=f"FAILED · {label}", border_style="red"))
        raise subprocess.CalledProcessError(code, script)
    seconds = time.perf_counter() - started
    console.print(f"  [green bold]✓ PASS[/]  {seconds:.1f} s")
    return seconds


def main():
    console = Console()
    console.print(Panel("Exact arithmetic · Python / SymPy–FLINT", title="ψ₃,₀ · recurrence certificate", border_style="cyan"))
    started = time.perf_counter()
    try:
        with (HERE / "python_runner.log").open("w") as log:
            recurrence = run_stage("run_certificate.py", "Integral identities and eight recurrence sectors", console, log)
            links = run_stage("check_inputs.py", "Coefficient-table links", console, log)
    except (subprocess.CalledProcessError, OSError) as exc:
        console.print(Text(f"Certificate failed: {exc}", style="bold red"))
        console.print("Details: python_runner.log and the individual checker logs.")
        return 1
    table = Table(header_style="bold", box=None, padding=(0, 2))
    table.add_column("Check")
    table.add_column("Result")
    table.add_column("Time", justify="right")
    table.add_row("Integral identities + 8/8 sectors", "[green]PASS[/]", f"{recurrence:.1f} s")
    table.add_row("Coefficient-table links", "[green]PASS[/]", f"{links:.1f} s")
    console.print(table)
    console.print(f"\n[bold green]✓ Certificate verified[/]  ·  {time.perf_counter()-started:.1f} s")
    console.print("[dim]Report: certificate.json · Details: python_runner.log[/]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
