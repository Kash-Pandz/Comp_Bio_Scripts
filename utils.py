from loguru import logger
import subprocess

def configure_logging(verbose: bool):
    """Configure loguru logging."""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(lambda msg: print(msg, end=""), level=level)


def run(command, description=""):
    """Run the shell command and stream stdout/stderr in real time."""
    logger.info(f"[RUN] {description}: {command}\n")
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        texts=True
    )
    for line in iter(process.stdout.readline, ''):
        if line:
            logger.info(line.strip())
    process.stdout.close()
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)
    return return_code