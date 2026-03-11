"""
FasTAN installer (adapted from satellome)
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

FASTAN_REPO = "https://github.com/ad3002/FASTAN.git"


def install_fastan(force: bool = False) -> bool:
    logger.info("Starting FasTAN installation...")

    bin_dir = get_bin_dir()
    fastan_path = bin_dir / 'fastan'

    if fastan_path.exists() and not force:
        logger.info(f"FasTAN already installed at {fastan_path}")
        if verify_installation('fastan'):
            return True
        logger.warning("Existing FasTAN binary failed verification, reinstalling...")

    deps_ok, error_msg = check_build_dependencies()
    if not deps_ok:
        logger.error(f"Build dependencies check failed:\n{error_msg}")
        return False

    with tempfile.TemporaryDirectory(prefix='fastan_build_') as tmp_dir:
        repo_dir = Path(tmp_dir) / 'FASTAN'

        try:
            logger.info(f"Cloning FasTAN repository...")
            result = subprocess.run(
                ['git', 'clone', FASTAN_REPO, str(repo_dir)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300
            )
            if result.returncode != 0:
                logger.error(f"Failed to clone FasTAN:\n{result.stderr.decode()}")
                return False

            # Patch Makefile for pthread
            makefile_path = repo_dir / 'Makefile'
            if makefile_path.exists():
                with open(makefile_path, 'r') as f:
                    content = f.read()
                if '-lpthread' not in content:
                    if 'LDFLAGS' in content:
                        content = content.replace('LDFLAGS =', 'LDFLAGS = -lpthread')
                    else:
                        content = 'LDFLAGS = -lpthread\n\n' + content
                    lines = content.split('\n')
                    modified = []
                    for line in lines:
                        if (('gcc' in line or 'cc' in line or '$(CC)' in line) and
                                '-o' in line and '$(LDFLAGS)' not in line and
                                not line.strip().startswith('#')):
                            line = line.replace(' -o ', ' $(LDFLAGS) -o ')
                        modified.append(line)
                    content = '\n'.join(modified)
                    with open(makefile_path, 'w') as f:
                        f.write(content)

            logger.info("Compiling FasTAN...")
            success, error_msg = run_make_with_fallback(repo_dir)
            if not success:
                logger.error(f"Failed to compile FasTAN:\n{error_msg}")
                return False

            # Find binary
            for search_dir in [repo_dir, repo_dir / 'bin']:
                if not search_dir.exists():
                    continue
                for name in ['FasTAN', 'fastan', 'FASTAN']:
                    candidate = search_dir / name
                    if candidate.exists() and os.access(candidate, os.X_OK):
                        shutil.copy2(candidate, fastan_path)
                        os.chmod(fastan_path, 0o755)
                        logger.info(f"FasTAN installed to {fastan_path}")
                        return verify_installation('fastan')

            logger.error(f"Could not find FasTAN binary after compilation")
            return False

        except subprocess.TimeoutExpired:
            logger.error("Installation timed out")
            return False


def uninstall_fastan() -> bool:
    fastan_path = get_bin_dir() / 'fastan'
    if not fastan_path.exists():
        return True
    fastan_path.unlink()
    logger.info("FasTAN uninstalled")
    return True
