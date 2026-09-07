"""List or selectively download members from a remote, uncompressed TAR.

The RADAR archive is multi-gigabyte, but its HTTP endpoint supports byte ranges.
This utility reads only TAR headers while listing and requests only the payload
range of selected members.  It therefore avoids downloading the whole archive.
"""

from __future__ import annotations

import argparse
import functools
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError
from urllib.request import Request, urlopen


RADAR_E35_URL = (
    "https://radar.kit.edu/radar-backend/archives/"
    "nWczysTWjmuzgFKm/versions/1/content"
)
BLOCK = 512
FETCH_CHUNK = 1024 * 1024


@dataclass(frozen=True)
class TarMember:
    name: str
    size: int
    data_offset: int
    typeflag: str


@functools.lru_cache(maxsize=4)
def _fetch_chunk(url: str, chunk_index: int) -> bytes:
    start = chunk_index * FETCH_CHUNK
    end = start + FETCH_CHUNK - 1
    request = Request(
        url,
        headers={
            "Range": f"bytes={start}-{end}",
            "User-Agent": "degali-selective-radar-reader/1.0",
        },
    )
    for attempt in range(8):
        try:
            with urlopen(request, timeout=60) as response:
                return response.read()
        except HTTPError as error:
            if error.code != 429 or attempt == 7:
                raise
            # The public archive throttles bursts but does not consistently
            # return Retry-After. Short repeated waits keep each pause bounded.
            time.sleep(10.0)
    raise AssertionError("unreachable")


def _range(url: str, start: int, end: int) -> bytes:
    pieces = []
    cursor = start
    while cursor <= end:
        chunk_index = cursor // FETCH_CHUNK
        chunk = _fetch_chunk(url, chunk_index)
        chunk_start = chunk_index * FETCH_CHUNK
        local_start = cursor - chunk_start
        local_end = min(len(chunk), end - chunk_start + 1)
        if local_start >= len(chunk) or local_end <= local_start:
            break
        pieces.append(chunk[local_start:local_end])
        cursor = chunk_start + local_end
    data = b"".join(pieces)
    expected = end - start + 1
    if len(data) != expected:
        raise OSError(
            f"range {start}-{end} returned {len(data)} bytes, expected {expected}"
        )
    return data


def _field(block: bytes, start: int, length: int) -> str:
    return block[start : start + length].split(b"\0", 1)[0].decode(
        "utf-8", errors="replace"
    )


def _octal(block: bytes, start: int, length: int) -> int:
    raw = block[start : start + length].strip(b"\0 ")
    return int(raw or b"0", 8)


def _pax_path(payload: bytes) -> str | None:
    position = 0
    while position < len(payload):
        space = payload.find(b" ", position)
        if space < 0:
            break
        record_length = int(payload[position:space])
        record = payload[space + 1 : position + record_length].rstrip(b"\n")
        key, separator, value = record.partition(b"=")
        if separator and key == b"path":
            return value.decode("utf-8", errors="replace")
        position += record_length
    return None


def iter_members(url: str):
    offset = 0
    pending_name: str | None = None
    while True:
        block = _range(url, offset, offset + BLOCK - 1)
        if block == b"\0" * BLOCK:
            return

        name = _field(block, 0, 100)
        prefix = _field(block, 345, 155)
        if prefix:
            name = f"{prefix}/{name}"
        size = _octal(block, 124, 12)
        typeflag = chr(block[156]) if block[156] else "0"
        data_offset = offset + BLOCK
        padded_size = ((size + BLOCK - 1) // BLOCK) * BLOCK

        if typeflag in {"L", "x"}:
            payload = _range(url, data_offset, data_offset + size - 1) if size else b""
            if typeflag == "L":
                pending_name = payload.rstrip(b"\0\n").decode(
                    "utf-8", errors="replace"
                )
            else:
                pending_name = _pax_path(payload) or pending_name
        else:
            yield TarMember(pending_name or name, size, data_offset, typeflag)
            pending_name = None

        offset = data_offset + padded_size


def safe_output_path(output_dir: Path, member_name: str) -> Path:
    leaf = PurePosixPath(member_name).name
    if not leaf or leaf in {".", ".."}:
        raise ValueError(f"unsafe TAR member name: {member_name!r}")
    return output_dir / leaf


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default=RADAR_E35_URL)
    parser.add_argument(
        "--match", help="case-insensitive regular expression for members to download"
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    pattern = re.compile(args.match, re.IGNORECASE) if args.match else None
    if pattern and args.output_dir is None:
        parser.error("--output-dir is required with --match")

    selected = 0
    for member in iter_members(args.url):
        print(f"{member.size:>12}  {member.typeflag}  {member.name}")
        if not pattern or member.typeflag not in {"0", "\0"}:
            continue
        if not pattern.search(member.name):
            continue
        output = safe_output_path(args.output_dir, member.name)
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = (
            _range(
                args.url,
                member.data_offset,
                member.data_offset + member.size - 1,
            )
            if member.size
            else b""
        )
        output.write_bytes(payload)
        print(f"downloaded -> {output}", file=sys.stderr)
        selected += 1

    if pattern and not selected:
        print("no matching regular file found", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
