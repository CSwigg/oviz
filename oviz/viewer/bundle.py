"""Binary bundle container used by the Oviz viewer.

A bundle is a JSON *manifest* plus a list of binary *blobs*. Blobs hold
columnar numeric arrays (positions over time, per-object attributes, volume
voxels, member-star catalogues) and opaque media (JPEG/PNG images). The
manifest references blobs by id and never contains bulk numeric data.

In the HTML file every blob becomes one ``<script type="application/
octet-stream">`` element holding base64 text. The browser decodes them with
native primitives (``Uint8Array.fromBase64`` or a ``data:`` fetch, then
``DecompressionStream("gzip")``), in priority order, so the scene can render
before large volumes finish arriving.

Encodings
---------
``gzip``          plain gzip of the little-endian array bytes.
``gzip+shuffle``  the array's bytes are byte-transposed (all first bytes of
                  every element, then all second bytes, ...) before gzip.
                  For smooth float data this typically shrinks the
                  compressed size by 20-40%; the browser un-shuffles.
``raw``           stored as-is (already-compressed media).
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import zlib
from dataclasses import dataclass, field
from typing import Any

import numpy as np

BUNDLE_FORMAT = "oviz-bundle/1"

_DTYPES = {
    "f32": np.float32,
    "f64": np.float64,
    "u8": np.uint8,
    "u16": np.uint16,
    "u32": np.uint32,
    "i8": np.int8,
    "i16": np.int16,
    "i32": np.int32,
}

# Lower numbers load first. The viewer starts rendering once every
# CRITICAL blob is decoded; the rest stream in afterwards.
CRITICAL = 0
NORMAL = 1
DEFERRED = 2


@dataclass
class Blob:
    id: str
    payload: bytes
    dtype: str
    shape: tuple[int, ...]
    codec: str
    priority: int = NORMAL
    mime: str | None = None
    raw_bytes: int = 0

    def descriptor(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "id": self.id,
            "dtype": self.dtype,
            "shape": list(self.shape),
            "codec": self.codec,
            "priority": self.priority,
            "bytes": len(self.payload),
            "rawBytes": self.raw_bytes,
        }
        if self.mime:
            d["mime"] = self.mime
        return d


def _gzip(data: bytes, level: int = 9) -> bytes:
    # mtime=0 keeps output byte-for-byte reproducible.
    return gzip.compress(data, compresslevel=level, mtime=0)


def byte_shuffle(data: bytes, itemsize: int) -> bytes:
    if itemsize <= 1 or not data:
        return data
    arr = np.frombuffer(data, dtype=np.uint8).reshape(-1, itemsize)
    return np.ascontiguousarray(arr.T).tobytes()


def byte_unshuffle(data: bytes, itemsize: int) -> bytes:
    if itemsize <= 1 or not data:
        return data
    arr = np.frombuffer(data, dtype=np.uint8).reshape(itemsize, -1)
    return np.ascontiguousarray(arr.T).tobytes()


@dataclass
class BundleBuilder:
    """Accumulates blobs, de-duplicating identical content."""

    blobs: list[Blob] = field(default_factory=list)
    _by_hash: dict[str, str] = field(default_factory=dict)
    _counter: int = 0
    compress_level: int = 9

    def _next_id(self, hint: str) -> str:
        self._counter += 1
        safe = "".join(c if c.isalnum() or c in "-_" else "-" for c in hint)[:40]
        return f"b{self._counter}-{safe}" if safe else f"b{self._counter}"

    def _register(self, key: str, make) -> str:
        existing = self._by_hash.get(key)
        if existing is not None:
            return existing
        blob = make()
        self.blobs.append(blob)
        self._by_hash[key] = blob.id
        return blob.id

    def array(
        self,
        values: Any,
        dtype: str,
        *,
        hint: str = "",
        priority: int = NORMAL,
        shuffle: bool | None = None,
    ) -> str:
        np_dtype = _DTYPES[dtype]
        arr = np.ascontiguousarray(np.asarray(values, dtype=np_dtype))
        data = arr.astype(np.dtype(np_dtype).newbyteorder("<"), copy=False).tobytes()
        itemsize = arr.dtype.itemsize
        if shuffle is None:
            shuffle = itemsize > 1 and arr.size >= 64
        key = hashlib.sha1(data + dtype.encode() + str(arr.shape).encode() + bytes([shuffle])).hexdigest()

        def make() -> Blob:
            body = byte_shuffle(data, itemsize) if shuffle else data
            return Blob(
                id=self._next_id(hint),
                payload=_gzip(body, self.compress_level),
                dtype=dtype,
                shape=tuple(int(s) for s in arr.shape),
                codec="gzip+shuffle" if shuffle else "gzip",
                priority=priority,
                raw_bytes=len(data),
            )

        return self._register(key, make)

    def media(self, data: bytes, mime: str, *, hint: str = "", priority: int = NORMAL) -> str:
        key = hashlib.sha1(data + mime.encode()).hexdigest()

        def make() -> Blob:
            return Blob(
                id=self._next_id(hint),
                payload=bytes(data),
                dtype="bytes",
                shape=(len(data),),
                codec="raw",
                priority=priority,
                mime=mime,
                raw_bytes=len(data),
            )

        return self._register(key, make)

    def json(self, obj: Any, *, hint: str = "", priority: int = NORMAL) -> str:
        data = json.dumps(obj, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")
        key = hashlib.sha1(data + b"json").hexdigest()

        def make() -> Blob:
            return Blob(
                id=self._next_id(hint),
                payload=_gzip(data, self.compress_level),
                dtype="json",
                shape=(len(data),),
                codec="gzip",
                priority=priority,
                raw_bytes=len(data),
            )

        return self._register(key, make)


@dataclass
class Bundle:
    manifest: dict[str, Any]
    blobs: list[Blob]

    def blob(self, blob_id: str) -> Blob:
        for b in self.blobs:
            if b.id == blob_id:
                return b
        raise KeyError(blob_id)

    def decode(self, blob_id: str) -> np.ndarray | bytes | Any:
        """Decode a blob back to numpy/bytes/JSON (used by tests and tools)."""

        b = self.blob(blob_id)
        if b.codec == "raw":
            return b.payload
        body = zlib.decompress(b.payload, 16 + zlib.MAX_WBITS)
        if b.dtype == "json":
            return json.loads(body.decode("utf-8"))
        np_dtype = np.dtype(_DTYPES[b.dtype]).newbyteorder("<")
        if b.codec == "gzip+shuffle":
            body = byte_unshuffle(body, np_dtype.itemsize)
        return np.frombuffer(body, dtype=np_dtype).reshape(b.shape)

    def size_report(self) -> dict[str, Any]:
        rows = sorted(self.blobs, key=lambda b: -len(b.payload))
        total = sum(len(b.payload) for b in self.blobs)
        return {
            "blobs": len(self.blobs),
            "payload_bytes": total,
            "base64_bytes": sum(4 * ((len(b.payload) + 2) // 3) for b in self.blobs),
            "largest": [(b.id, len(b.payload), b.raw_bytes) for b in rows[:12]],
        }


def encode_base64_lines(data: bytes, width: int = 1 << 16) -> str:
    text = base64.b64encode(data).decode("ascii")
    if len(text) <= width:
        return text
    return "\n".join(text[i:i + width] for i in range(0, len(text), width))
