"""Sequence helpers.

Conventions (fixed on day 1, never change them silently):
  * base codes: A=0 C=1 G=2 T/U=3, anything else = 4 (never matches anything)
  * all internal coordinates are 0-based, half-open [start, end)
  * output files use 1-based inclusive coordinates and say so in the header
"""
import gzip
import hashlib

import numpy as np

_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _c in enumerate("ACGT"):
    _LUT[ord(_c)] = _i
    _LUT[ord(_c.lower())] = _i
_LUT[ord("U")] = 3
_LUT[ord("u")] = 3

SEP = 4  # separator / ambiguity code


def encode(s):
    if isinstance(s, str):
        s = s.encode("ascii")
    return _LUT[np.frombuffer(s, dtype=np.uint8)]


def decode(arr):
    return "".join("ACGTN"[min(int(x), 4)] for x in arr)


def revcomp(arr):
    out = arr[::-1].copy()
    ok = out < 4
    out[ok] = 3 - out[ok]
    return out


def _open(path):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")


def read_fasta(path):
    """Yield (name, raw_bytes) for every record."""
    name, parts = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(parts)
                name, parts = line[1:].split()[0].decode(), []
            elif line:
                parts.append(line)
    if name is not None:
        yield name, b"".join(parts)


class Reference:
    """A whole FASTA as one uint8 array, records separated by SEP so no match
    can span two transcripts. `locate` maps a global position back to a record."""

    def __init__(self, seq, starts, names):
        self.seq, self.starts, self.names = seq, starts, names

    @classmethod
    def from_fasta(cls, path):
        arrays, starts, names, pos = [], [], [], 0
        gap = np.array([SEP], dtype=np.uint8)
        for name, raw in read_fasta(path):
            a = encode(raw)
            starts.append(pos)
            names.append(name)
            arrays += [a, gap]
            pos += a.size + 1
        if not arrays:
            raise ValueError(f"no FASTA records in {path}")
        return cls(np.concatenate(arrays), np.array(starts, dtype=np.int64), names)

    @classmethod
    def from_string(cls, s, name="seq"):
        a = encode(s)
        return cls(a, np.array([0], dtype=np.int64), [name])

    def locate(self, pos):
        i = int(np.searchsorted(self.starts, pos, side="right")) - 1
        return self.names[i], int(pos - self.starts[i])


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()
