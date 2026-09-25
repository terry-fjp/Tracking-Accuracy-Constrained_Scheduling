"""Verify all released scene fields and exact numerical provenance."""
import hashlib
import json
from .paths import DATA
from .scenes import load_scene

def main():
    manifest=json.loads((DATA/"scenes/manifest.json").read_text())
    for item in manifest:
        model=load_scene(DATA/"scenes"/f"{item['scene_id']}.json")
        digest=hashlib.sha256(json.dumps(model,sort_keys=True,separators=(",",":")).encode()).hexdigest()
        if digest!=item["canonical_model_sha256"]:
            raise ValueError(f"Scene {item['scene_id']}: reference numerical content changed")
        if len(model["targets"])!=item["targets"] or len(model["observers"])!=item["observers"]:
            raise ValueError("Scene count mismatch")
    print(json.dumps({"status":"PASS","scenes":len(manifest),"targets":sum(x["targets"] for x in manifest),"numerical_payload_preserved":True}))
