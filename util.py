from subprocess import Popen, call
import traceback, re
import logging

MODULE_LOGGER = logging.getLogger('prepare_data_construct_tree.util')

REGEX_VOID = re.compile(r"(\(|\)|:|,|}|{|'|/|]|\[|\\)")
REGEX_UNDERSCORE = re.compile(r"( |\|)")
REGEX_VOID_SUBST = ""
REGEX_UNDERSCORE_SUBST = "_"

def prepareNames(line):
	line = REGEX_VOID.sub(REGEX_VOID_SUBST, line)
	return REGEX_UNDERSCORE.sub(REGEX_UNDERSCORE_SUBST, line)
	
def runSubProcess(command):
	try:
		call(command, shell=True)
	except OSError, osError:
		MODULE_LOGGER.error("osError " + osError)
		MODULE_LOGGER.error(traceback.print_exc())
	except Exception, e:
		MODULE_LOGGER.error(traceback.print_exc())
		
def runSubProcessWithCheck(command, processName):
	try:
		proc = Popen(command, shell=True)
		status = proc.poll()
		while status == None:
			MODULE_LOGGER.info("Still runnig " + processName)
			status = proc.poll()
		MODULE_LOGGER.info("Seems like finished " + processName + " " + str(status))
	except OSError, osError:
		MODULE_LOGGER.error("osError " + osError)
		MODULE_LOGGER.error(traceback.print_exc())
	except Exception, e:
		MODULE_LOGGER.error(traceback.print_exc())
