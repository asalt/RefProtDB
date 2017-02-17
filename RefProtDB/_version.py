import os
from datetime import date

BASE_DIR = os.path.dirname(os.path.realpath(__file__))

def _calculate_latest_edit(files, path=BASE_DIR):
    return max([os.stat(os.path.join(path, f)).st_mtime for f in files])

_package_files = ('_version.py', 'cli.py', 'containers.py', 'download_utils.py', 'utils.py')
__version__ = '0.1.0'
__copyright__ = date.fromtimestamp(_calculate_latest_edit(_package_files)).__str__()
