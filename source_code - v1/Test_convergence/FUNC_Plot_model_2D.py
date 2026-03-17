import importlib.util
from pathlib import Path

_MODULE_PATH = Path(__file__).resolve().parent.parent / "FUNC_Plot_model_2D.py"
_SPEC = importlib.util.spec_from_file_location("base_FUNC_Plot_model_2D", _MODULE_PATH)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

FUNC_Plot_model_2D = _MODULE.FUNC_Plot_model_2D
