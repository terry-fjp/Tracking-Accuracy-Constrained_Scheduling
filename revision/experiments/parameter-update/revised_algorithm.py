"""Revised warning policy; original frozen algorithm and results stay unchanged."""
from pathlib import Path
import json,copy
CONFIG=Path(__file__).with_name('algorithm-config.json')
def revised_scene(scene):
    out=copy.deepcopy(scene)
    out['config']['warning']=json.loads(CONFIG.read_text())['warning']
    return out
