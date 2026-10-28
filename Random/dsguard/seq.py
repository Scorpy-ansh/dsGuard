import gzip

import hashlib



import numpy as np



_LUT = np.full(256, 4, dtype=np.uint8)

for _i, _c in enumerate("ACGT"):

    _LUT[ord(_c)] = _i

    _LUT[ord(_c.lower())] = _i

_LUT[ord("U")] = 3

_LUT[ord("u")] = 3



SEP = 4





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







def _iter_records(path):

    name, description, parts = None, "", []



    with _open(path) as fh:

        for line_number, line in enumerate(fh, start=1):

            line = line.strip()



            if not line:

                continue



            if line.startswith(b">"):

                if name is not None:

                    yield name, description, b"".join(parts)



                header = line[1:].split(None, 1)

                if not header:

                    raise ValueError(

                        f"Invalid FASTA header at line {line_number}: missing sequence name"

                    )



                name = header[0].decode("ascii")

                description = (

                    header[1].decode("ascii", "replace").strip()

                    if len(header) > 1 else ""

                )

                parts = []



            else:

                if name is None:

                    raise ValueError(

                        f"Invalid FASTA at line {line_number}: sequence data before first header"

                    )



                parts.append(line)



    if name is not None:

        yield name, description, b"".join(parts)





def read_fasta(path):

    for name, _description, raw in _iter_records(path):

        yield name, raw





class Reference:



    def __init__(self, seq, starts, names, descriptions=None):

        self.seq, self.starts, self.names = seq, starts, names

        self.descriptions = (

            list(descriptions) if descriptions is not None else [""] * len(names)

        )



    @classmethod

    def from_fasta(cls, path):

        arrays, starts, names, descriptions, pos = [], [], [], [], 0

        gap = np.array([SEP], dtype=np.uint8)

        for name, description, raw in _iter_records(path):

            a = encode(raw)

            starts.append(pos)

            names.append(name)

            descriptions.append(description)

            arrays += [a, gap]

            pos += a.size + 1

        if not arrays:

            raise ValueError(f"no FASTA records in {path}")

        return cls(np.concatenate(arrays), np.array(starts, dtype=np.int64),

                   names, descriptions)



    @classmethod

    def from_string(cls, s, name="seq"):

        a = encode(s)

        return cls(a, np.array([0], dtype=np.int64), [name], [""])



    @property

    def records_without_annotation(self):

        return sum(1 for description in self.descriptions if not description)



    def locate(self, pos):

        i = int(np.searchsorted(self.starts, pos, side="right")) - 1

        return self.names[i], int(pos - self.starts[i])





def sha256_file(path):

    h = hashlib.sha256()

    with open(path, "rb") as fh:

        for block in iter(lambda: fh.read(1 << 20), b""):

            h.update(block)

    return h.hexdigest()
