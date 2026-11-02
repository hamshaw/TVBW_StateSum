import sys
from .sageRegina.test import runTests

if __name__ == '__main__':
    success = runTests()
    sys.exit(not success)
