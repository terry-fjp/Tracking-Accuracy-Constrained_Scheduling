"""Verify and prepare the frozen revision experiments using only the standard library."""
import argparse
import hashlib
import json
import shutil
import stat
import zipfile
from pathlib import Path, PurePosixPath

REPOSITORY = Path(__file__).resolve().parents[1]
ROOT = REPOSITORY / 'revision/experiments'
FINGERPRINTS = {
    'protocol.json': '897fc155a0a40c679acec45d26910c621b77295215516c8d012c1cbfc0a2d975',
    'parameter-update/protocol.json': 'a7c65e81362b54812bc5411f24e5c875691ab46993f07bc152a678358e85590d',
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe_target(base, name):
    path = PurePosixPath(name)
    require(bool(name) and not path.is_absolute() and '..' not in path.parts
            and '\\' not in name and ':' not in name, f'Unsafe archive path: {name}')
    target = base.joinpath(*path.parts)
    for node in [target, *target.parents]:
        require(not node.is_symlink(), f'Symlink in destination: {node}')
        if node == REPOSITORY:
            break
    require(target.resolve().is_relative_to(base.resolve()), f'Path outside archive: {name}')
    return target


def prepare(check_only=False):
    protocols = {}
    for relative, expected in FINGERPRINTS.items():
        p = json.loads((ROOT / relative).read_text())
        fingerprint = p.pop('fingerprint')
        require(fingerprint == expected == sha(json.dumps(p, sort_keys=True).encode()),
                f'Frozen protocol changed: {relative}')
        protocols[relative] = p
    parent = protocols['protocol.json']
    validation = protocols['parameter-update/protocol.json']
    for relative, expected in parent['source_hashes'].items():
        require(sha((ROOT / relative).read_bytes()) == expected, f'Frozen source changed: {relative}')
    for relative, expected in validation['code_hashes'].items():
        require(sha((ROOT/'parameter-update'/relative).read_bytes()) == expected,
                f'Validation code changed: {relative}')
    require(validation['parent_physical_protocol'] == FINGERPRINTS['protocol.json'],
            'Validation parent protocol mismatch')
    selection = ROOT/'reference-results/summary.json'
    require(sha(selection.read_bytes()) == validation['selection_data_sha256'],
            'Frozen parameter-selection record changed')
    active_selection = safe_target(ROOT, 'results/summary.json')
    if not check_only and active_selection.exists():
        require(sha(active_selection.read_bytes()) == validation['selection_data_sha256'],
                'Existing results/summary.json differs from the frozen selection record. '
                'Use a fresh checkout for validation; no output has been overwritten.')
    manifest = json.loads((ROOT/'provenance/reference-archive-hashes.json').read_text())
    base = ROOT/'reference-archive'
    with zipfile.ZipFile(ROOT/'reference-archive.zip') as archive:
        entries = archive.infolist()
        require(len(entries) == len(manifest) == 1339, 'Unexpected archive file count')
        require(len({m.filename for m in entries}) == len(entries), 'Duplicate archive entry')
        require({m.filename for m in entries} == set(manifest), 'Archive/manifest paths differ')
        for member in entries:
            target = safe_target(base, member.filename)
            require(not member.is_dir() and not stat.S_ISLNK(member.external_attr >> 16),
                    f'Non-file archive entry: {member.filename}')
            data = archive.read(member)
            require(sha(data) == manifest[member.filename], f'Archive hash mismatch: {member.filename}')
            if target.exists():
                require(target.is_file() and sha(target.read_bytes()) == manifest[member.filename],
                        f'Existing baseline differs: {target}')
        if not check_only:
            for member in entries:
                target = safe_target(base, member.filename)
                if not target.exists():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with target.open('xb') as stream:
                        stream.write(archive.read(member))
    if not check_only:
        for relative in ['checks', 'results', 'parameter-update/results', '../paper/figures', '../audit']:
            (ROOT/relative).mkdir(parents=True, exist_ok=True)
        if not active_selection.exists():
            shutil.copyfile(selection, active_selection)
    print(json.dumps({
        'status': 'PASS', 'mode': 'check-only' if check_only else 'prepared',
        'frozen_source_files': len(parent['source_hashes']),
        'validation_code_files': len(validation['code_hashes']),
        'baseline_archive_files': len(manifest),
        'selection_sha256': validation['selection_data_sha256'],
        'scope': 'Input integrity and workspace preparation; no scheduling or Monte Carlo experiment executed.'
    }, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true', help='Verify inputs without extracting or creating outputs.')
    args = parser.parse_args()
    try:
        prepare(args.check_only)
    except (ValueError, OSError, KeyError, zipfile.BadZipFile) as error:
        parser.exit(1, f'Preparation failed: {error}\n')
