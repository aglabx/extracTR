"""
Base utilities for installers (adapted from satellome)
"""

import os
import shutil
import platform
import logging
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


def detect_platform() -> Tuple[str, str]:
    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == 'linux':
        platform_name = 'linux'
    elif system == 'darwin':
        platform_name = 'darwin'
    else:
        platform_name = 'unknown'

    if machine in ['x86_64', 'amd64']:
        arch = 'x86_64'
    elif machine in ['arm64', 'aarch64']:
        arch = 'arm64'
    else:
        arch = 'unknown'

    return platform_name, arch


def get_bin_dir() -> Path:
    """Get or create the extracTR binary directory."""
    try:
        import extractr
        package_dir = Path(extractr.__file__).parent
        bin_dir = package_dir / 'bin'
        bin_dir.mkdir(parents=True, exist_ok=True)
        test_file = bin_dir / '.write_test'
        try:
            test_file.touch()
            test_file.unlink()
            return bin_dir
        except (PermissionError, OSError):
            pass
    except Exception:
        pass

    bin_dir = Path.home() / '.extractr' / 'bin'
    bin_dir.mkdir(parents=True, exist_ok=True)
    return bin_dir


def check_binary_exists(binary_name: str) -> Optional[str]:
    """Check if a binary exists in extracTR bin directory or system PATH."""
    extractr_bin = get_bin_dir() / binary_name
    if extractr_bin.exists() and os.access(extractr_bin, os.X_OK):
        return str(extractr_bin)
    system_path = shutil.which(binary_name)
    if system_path:
        return system_path
    return None


def verify_installation(binary_name: str) -> bool:
    binary_path = check_binary_exists(binary_name)
    if not binary_path:
        logger.error(f"{binary_name} not found")
        return False
    logger.info(f"Found {binary_name} at: {binary_path}")
    return True


def find_system_gcc() -> Optional[str]:
    for gcc_path in ['/usr/bin/gcc', '/usr/local/bin/gcc']:
        if os.path.exists(gcc_path):
            return gcc_path
    return None


def run_make_with_fallback(cwd: Path, timeout: int = 300) -> Tuple[bool, str]:
    import subprocess

    compile_env = os.environ.copy()

    result = subprocess.run(
        ['make'], cwd=cwd,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=timeout, env=compile_env
    )

    if result.returncode == 0:
        return True, ""

    stderr_output = result.stderr.decode()

    if 'cannot find -lz' in stderr_output and 'conda' in os.environ.get('PATH', ''):
        logger.warning("Compilation failed with conda gcc (cannot find -lz), retrying with system gcc...")
        subprocess.run(['make', 'clean'], cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        system_gcc = find_system_gcc()
        if system_gcc:
            compile_env['CC'] = system_gcc
            result = subprocess.run(
                ['make'], cwd=cwd,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=timeout, env=compile_env
            )
            if result.returncode == 0:
                return True, ""
            return False, f"Failed with system gcc:\n{result.stderr.decode()}"
        return False, f"System gcc not found. Original error:\n{stderr_output}"

    return False, f"Compilation failed:\n{stderr_output}"


def check_build_dependencies() -> Tuple[bool, str]:
    missing = []
    if not shutil.which('git'):
        missing.append('git')
    if not shutil.which('make'):
        missing.append('make')
    if not (shutil.which('gcc') or shutil.which('clang') or shutil.which('cc')):
        missing.append('gcc or clang')

    if missing:
        error_msg = f"Missing build tools: {', '.join(missing)}\n"
        system_name, _ = detect_platform()
        if system_name == 'darwin':
            error_msg += "On macOS: xcode-select --install"
        elif system_name == 'linux':
            error_msg += "On Ubuntu/Debian: sudo apt-get install build-essential git"
        return False, error_msg

    return True, ""
