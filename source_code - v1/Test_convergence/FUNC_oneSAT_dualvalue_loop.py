import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from FUNC_oneSAT_dualvalue_loop import FUNC_oneSAT_dualvalue_loop  # noqa: E402,F401
