"""Lanzador local y deliberadamente pequeño para la aplicación Marimo."""

from __future__ import annotations

import argparse
import socket
import subprocess
import sys
from pathlib import Path
from typing import Sequence

DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 2718
MINIMUM_PORT = 1024
MAXIMUM_PORT = 65535
APPLICATION = Path(__file__).resolve().parents[1] / "ui" / "gps_app.py"


def local_port(value: str) -> int:
    """Parse a non-privileged TCP port for the loopback listener."""
    try:
        port = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("el puerto debe ser un entero") from error
    if not MINIMUM_PORT <= port <= MAXIMUM_PORT:
        raise argparse.ArgumentTypeError(
            f"el puerto debe pertenecer a {MINIMUM_PORT}..={MAXIMUM_PORT}"
        )
    return port


def port_is_available(port: int) -> bool:
    """Check whether the requested loopback port can be bound now."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
        try:
            listener.bind((DEFAULT_HOST, port))
        except OSError:
            return False
    return True


def marimo_command(port: int, no_open: bool) -> list[str]:
    """Build the argument-vector invocation without using a shell."""
    command = [
        sys.executable,
        "-m",
        "marimo",
        "run",
        str(APPLICATION),
        "--host",
        DEFAULT_HOST,
        "--port",
        str(port),
    ]
    if no_open:
        command.append("--headless")
    return command


def parse_args(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inicia la interfaz local de Genome Population Simulator."
    )
    parser.add_argument(
        "--port",
        type=local_port,
        default=DEFAULT_PORT,
        help=f"puerto local no privilegiado (por defecto: {DEFAULT_PORT})",
    )
    parser.add_argument(
        "--no-open",
        action="store_true",
        help="no abrir el navegador automáticamente",
    )
    return parser.parse_args(arguments)


def run_marimo(command: Sequence[str]) -> int:
    """Run Marimo and turn an interactive interruption into a clean shutdown."""
    process = subprocess.Popen(command)
    try:
        return process.wait()
    except KeyboardInterrupt:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        return 130


def main(arguments: Sequence[str] | None = None) -> int:
    options = parse_args(arguments)
    if not port_is_available(options.port):
        print(
            f"el puerto local {options.port} ya está ocupado; elija otro con --port",
            file=sys.stderr,
        )
        return 2

    url = f"http://{DEFAULT_HOST}:{options.port}"
    print(f"Interfaz local: {url}", flush=True)
    return run_marimo(marimo_command(options.port, options.no_open))


if __name__ == "__main__":
    raise SystemExit(main())
