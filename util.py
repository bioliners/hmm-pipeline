from subprocess import Popen, run, CalledProcessError
import re
import logging, time

MODULE_LOGGER = logging.getLogger('psi_blast_getmax.util')

REGEX_VOID = re.compile(r"(\(|\)|:|,|}|{|'|/|]|\[|\\)")
REGEX_UNDERSCORE = re.compile(r"( |\|)")
REGEX_VOID_SUBST = ""
REGEX_UNDERSCORE_SUBST = "_"

def prepareNames(line):
	line = REGEX_VOID.sub(REGEX_VOID_SUBST, line)
	return REGEX_UNDERSCORE.sub(REGEX_UNDERSCORE_SUBST, line)
	
def runSubProcess(command):
    try:
        run(command, shell=True, check=True)
    except CalledProcessError:
        MODULE_LOGGER.exception("Subprocess failed")
    except Exception:
        MODULE_LOGGER.exception("Unexpected error")
		
def runSubProcessWithCheck(command, processName, poll_interval=1.0):
    try:
        proc = Popen(command, shell=True)

        while True:
            status = proc.poll()
            if status is not None:
                break
            MODULE_LOGGER.info("Still running %s", processName)
            time.sleep(poll_interval)

        MODULE_LOGGER.info("Seems like finished %s %s", processName, status)
        return status

    except OSError:
        MODULE_LOGGER.exception("OSError while running process %s", processName)
        raise

    except Exception:
        MODULE_LOGGER.exception("Unexpected error while running process %s", processName)
        raise
