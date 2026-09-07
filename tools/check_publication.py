"""Fail fast when the proposed GitHub snapshot is unsafe or too large."""

from __future__ import annotations

import pathlib
import re
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
MAX_FILE = 50 * 1024 * 1024
MAX_TOTAL = 200 * 1024 * 1024
FORBIDDEN_SUFFIXES = {".xls", ".xlsx", ".dbf", ".pdf", ".pem", ".p12", ".pfx"}
FORBIDDEN_NAMES = {".env", "id_rsa", "id_ed25519", "credentials.json"}
SECRET_PATTERNS = {
    "private key": re.compile(rb"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
    "GitHub token": re.compile(rb"(?:ghp|github_pat)_[A-Za-z0-9_]{20,}"),
    "OpenAI-style key": re.compile(rb"\bsk-[A-Za-z0-9_-]{20,}"),
    "AWS access key": re.compile(rb"\bAKIA[0-9A-Z]{16}\b"),
}
ESSENTIAL = {
    "README.md",
    "LICENSE",
    "CITATION.cff",
    "pyproject.toml",
    ".github/workflows/tests.yml",
}


def publication_files() -> list[pathlib.Path]:
    command = [
        "git",
        "ls-files",
        "--cached",
        "--others",
        "--exclude-standard",
        "-z",
    ]
    result = subprocess.run(command, cwd=ROOT, check=False, capture_output=True)
    if result.returncode:
        raise RuntimeError("initialize the repository with git init before publication checks")
    names = [name for name in result.stdout.decode("utf-8").split("\0") if name]
    return [ROOT / name for name in names if (ROOT / name).is_file()]


def cff_version() -> str:
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r'^version:\s*["\']([^"\']+)', text, re.MULTILINE)
    if not match:
        raise ValueError("CITATION.cff has no version")
    return match.group(1)


def project_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    project = text.split("[project]", 1)[1]
    match = re.search(r'^version\s*=\s*["\']([^"\']+)', project, re.MULTILINE)
    if not match:
        raise ValueError("pyproject.toml [project] has no version")
    return match.group(1)


def main() -> int:
    paths = publication_files()
    relative = {path.relative_to(ROOT).as_posix() for path in paths}
    errors: list[str] = []
    warnings: list[str] = []

    missing = sorted(ESSENTIAL - relative)
    if missing:
        errors.append(f"essential public files missing: {missing}")

    total = sum(path.stat().st_size for path in paths)
    if total > MAX_TOTAL:
        errors.append(f"candidate snapshot is {total / 1024**2:.1f} MiB; limit is 200 MiB")

    for path in paths:
        rel = path.relative_to(ROOT).as_posix()
        size = path.stat().st_size
        if size > MAX_FILE:
            errors.append(f"file exceeds 50 MiB: {rel} ({size / 1024**2:.1f} MiB)")
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f"third-party/raw or credential file selected: {rel}")
        if path.name.lower() in FORBIDDEN_NAMES:
            errors.append(f"credential filename selected: {rel}")
        if size <= 5 * 1024 * 1024:
            data = path.read_bytes()
            for label, pattern in SECRET_PATTERNS.items():
                if pattern.search(data):
                    errors.append(f"possible {label} in {rel}")

    package_version = project_version()
    citation_version = cff_version()
    if package_version != citation_version:
        errors.append(
            f"version mismatch: pyproject={package_version}, CITATION={citation_version}"
        )

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    if "REPLACE-BEFORE-RELEASE" in citation:
        warnings.append("replace the CITATION.cff author placeholder before the first tag")
    if "github.com/lyullee/degali" not in citation:
        errors.append("CITATION.cff does not use the confirmed GitHub repository URL")

    if errors:
        print("PUBLICATION CHECK FAILED")
        for error in errors:
            print(f"- {error}")
        return 1

    print(
        f"PUBLICATION CHECK PASSED: {len(paths)} files, "
        f"{total / 1024**2:.1f} MiB, largest below 50 MiB"
    )
    for warning in warnings:
        print(f"MANUAL BEFORE TAG: {warning}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
