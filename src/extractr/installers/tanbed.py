"""
tanbed (from alntools) installer (adapted from satellome)
"""

import os
import shutil
import subprocess
import tempfile
import logging
from pathlib import Path

from .base import (
    get_bin_dir,
    check_build_dependencies,
    verify_installation,
    run_make_with_fallback
)

logger = logging.getLogger(__name__)

ALNTOOLS_REPO = "https://github.com/richarddurbin/alntools.git"


def install_tanbed(force: bool = False) -> bool:
    logger.info("Starting tanbed installation...")

    bin_dir = get_bin_dir()
    tanbed_path = bin_dir / 'tanbed'

    if tanbed_path.exists() and not force:
        logger.info(f"tanbed already installed at {tanbed_path}")
        if verify_installation('tanbed'):
            return True
        logger.warning("Existing tanbed binary failed verification, reinstalling...")

    deps_ok, error_msg = check_build_dependencies()
    if not deps_ok:
        logger.error(f"Build dependencies check failed:\n{error_msg}")
        return False

    with tempfile.TemporaryDirectory(prefix='tanbed_build_') as tmp_dir:
        repo_dir = Path(tmp_dir) / 'alntools'

        try:
            logger.info(f"Cloning alntools repository...")
            result = subprocess.run(
                ['git', 'clone', ALNTOOLS_REPO, str(repo_dir)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300
            )
            if result.returncode != 0:
                logger.error(f"Failed to clone alntools:\n{result.stderr.decode()}")
                return False

            logger.info("Compiling alntools...")
            success, error_msg = run_make_with_fallback(repo_dir)
            if not success:
                logger.error(f"Failed to compile alntools:\n{error_msg}")
                return False

            # Find tanbed binary
            for search_dir in [repo_dir, repo_dir / 'bin']:
                if not search_dir.exists():
                    continue
                candidate = search_dir / 'tanbed'
                if candidate.exists() and os.access(candidate, os.X_OK):
                    shutil.copy2(candidate, tanbed_path)
                    os.chmod(tanbed_path, 0o755)
                    logger.info(f"tanbed installed to {tanbed_path}")
                    return verify_installation('tanbed')

            logger.error(f"Could not find tanbed binary after compilation")
            return False

        except subprocess.TimeoutExpired:
            logger.error("Installation timed out")
            return False


def uninstall_tanbed() -> bool:
    tanbed_path = get_bin_dir() / 'tanbed'
    if not tanbed_path.exists():
        return True
    tanbed_path.unlink()
    logger.info("tanbed uninstalled")
    return True
