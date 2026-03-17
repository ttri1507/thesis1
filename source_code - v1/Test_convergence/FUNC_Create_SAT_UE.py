import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from FUNC_Create_SAT_UE import FUNC_Create_SAT_UE  # noqa: E402,F401
