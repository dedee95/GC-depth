#!/usr/bin/env python3
# Author: Dede Kurniawan
# Email: dedekurniawan@genomics.cn or dedearkun2710@gmail.com

"""
Compute and visualize GC content versus sequencing depth per genomic window.

Usage: gc_depth <fasta> <pandepth_output> [options]

Positional arguments:
  fasta                Genome FASTA file (gzipped is also fine)
  pandepth             Pandepth windowed depth file (.win.stat.gz)

Options:
  -h, --help           Show this help message and exit
  -w, --window WINDOW  Window size, must match pandepth -w value (default: 1000)
  -o, --output OUTPUT  Output plot file (.png or .pdf, default: gc-depth.png)
  --log-depth          Use logarithmic scale for the depth axis
  --plot-only TSV      Skip processing, re-plot from an existing combined TSV (from --output-data)
  --output-data FILE   Save merged GC and depth data to this TSV file (can be reused with --plot-only)
  --version            Show the installed version and exit
"""

import argparse
from array import array
import gzip
import math
import os
from pathlib import Path
import re
import sqlite3
import sys
import tempfile
from typing import BinaryIO, Callable, Dict, Optional, TextIO, Tuple

import matplotlib as mpl

# Use a non-interactive backend so the CLI works on servers and HPC systems.
mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.ndimage import gaussian_filter

from . import __version__


PAIR_BUFFER_FLOATS = 131_072
PLOT_CHUNK_SIZE = 250_000
DEPTH_RECORD_BUFFER = 65_536
FLOAT32_MAX = float(np.finfo(np.float32).max)
GC_LOOKUP = np.zeros(256, dtype=np.uint8)
GC_LOOKUP[ord("G")] = 1
GC_LOOKUP[ord("C")] = 1
INVALID_SEQUENCE_BYTES = re.compile(rb"[\x00-\x20\x7f>]")


class GCDepthError(Exception):
    """Raised for clear, user-facing input or processing errors."""


class DocstringHelpParser(argparse.ArgumentParser):
    """ArgumentParser that prints the module docstring verbatim for --help."""

    def format_help(self) -> str:
        doc = (__doc__ or "").strip()
        return f"{doc}\n" if doc else super().format_help()


def _display_path(path: os.PathLike) -> str:
    """Return a readable path without requiring the file to exist."""
    return str(Path(path).expanduser())


def open_file(path: os.PathLike) -> TextIO:
    """Open a UTF-8 plain-text or gzip-compressed file using magic-byte detection."""
    file_path = Path(path).expanduser()
    try:
        if not file_path.exists():
            raise GCDepthError(f"File not found: {_display_path(file_path)}")
        if not file_path.is_file():
            raise GCDepthError(f"Not a regular file: {_display_path(file_path)}")
        with file_path.open("rb") as raw:
            magic = raw.read(2)
    except PermissionError as exc:
        raise GCDepthError(f"Permission denied while opening: {_display_path(file_path)}") from exc
    except OSError as exc:
        raise GCDepthError(f"Cannot open {_display_path(file_path)}: {exc}") from exc

    try:
        if magic == b"\x1f\x8b":
            return gzip.open(file_path, "rt", encoding="utf-8", errors="strict", newline="")
        return file_path.open("rt", encoding="utf-8", errors="strict", newline="")
    except OSError as exc:
        raise GCDepthError(f"Cannot read {_display_path(file_path)}: {exc}") from exc


def open_binary_file(path: os.PathLike) -> BinaryIO:
    """Open a plain or gzip-compressed file in binary mode using magic-byte detection."""
    file_path = Path(path).expanduser()
    try:
        if not file_path.exists():
            raise GCDepthError(f"File not found: {_display_path(file_path)}")
        if not file_path.is_file():
            raise GCDepthError(f"Not a regular file: {_display_path(file_path)}")
        with file_path.open("rb") as raw:
            magic = raw.read(2)
        if magic == b"\x1f\x8b":
            return gzip.open(file_path, "rb")
        return file_path.open("rb")
    except PermissionError as exc:
        raise GCDepthError(f"Permission denied while opening: {_display_path(file_path)}") from exc
    except OSError as exc:
        raise GCDepthError(f"Cannot open {_display_path(file_path)}: {exc}") from exc


def _stream_fasta_events(
    path: os.PathLike,
    on_header: Callable[[str, int], None],
    on_sequence: Callable[[bytes, int], None],
    chunk_size: int = 1_048_576,
) -> bool:
    """Parse FASTA in fixed-size binary chunks, including unwrapped sequences."""
    line_number = 1
    mode: Optional[str] = None
    header_buffer = bytearray()
    saw_nonblank_content = False
    skip_lf_after_cr = False

    def finish_header() -> None:
        try:
            header_text = header_buffer.decode("utf-8", errors="strict")
        except UnicodeDecodeError as exc:
            raise GCDepthError(
                f"Cannot decode FASTA header at line {line_number} as UTF-8."
            ) from exc
        on_header(header_text, line_number)
        header_buffer.clear()

    try:
        with open_binary_file(path) as handle:
            while True:
                chunk = handle.read(chunk_size)
                if not chunk:
                    break
                position = 0
                chunk_length = len(chunk)

                if skip_lf_after_cr:
                    if chunk.startswith(b"\n"):
                        position = 1
                    skip_lf_after_cr = False

                while position < chunk_length:
                    if mode is None:
                        current_byte = chunk[position]
                        if current_byte in (10, 13):
                            position += 1
                            if current_byte == 13:
                                if position < chunk_length and chunk[position] == 10:
                                    position += 1
                                elif position == chunk_length:
                                    skip_lf_after_cr = True
                            line_number += 1
                            continue

                        saw_nonblank_content = True
                        if current_byte == 62:  # '>'
                            mode = "header"
                            header_buffer.clear()
                            position += 1
                        else:
                            mode = "sequence"

                    lf_position = chunk.find(b"\n", position)
                    cr_position = chunk.find(b"\r", position)
                    candidates = [value for value in (lf_position, cr_position) if value >= 0]
                    line_end = min(candidates) if candidates else chunk_length
                    segment = chunk[position:line_end]

                    if mode == "header":
                        header_buffer.extend(segment)
                    elif segment:
                        on_sequence(segment, line_number)

                    if line_end == chunk_length:
                        position = chunk_length
                        continue

                    line_break = chunk[line_end]
                    if mode == "header":
                        finish_header()
                    mode = None
                    line_number += 1
                    position = line_end + 1

                    if line_break == 13:
                        if position < chunk_length and chunk[position] == 10:
                            position += 1
                        elif position == chunk_length:
                            skip_lf_after_cr = True

            if mode == "header":
                finish_header()

    except (gzip.BadGzipFile, EOFError) as exc:
        raise GCDepthError(f"Cannot decompress FASTA file: {exc}") from exc

    return saw_nonblank_content


def _prepare_output_path(path: os.PathLike, label: str) -> Path:
    """Validate an output location before expensive processing begins."""
    output_path = Path(path).expanduser()
    parent = output_path.parent if str(output_path.parent) else Path(".")
    if not parent.exists():
        raise GCDepthError(f"{label} directory does not exist: {_display_path(parent)}")
    if not parent.is_dir():
        raise GCDepthError(f"{label} parent is not a directory: {_display_path(parent)}")
    if output_path.exists() and output_path.is_dir():
        raise GCDepthError(f"{label} path is a directory: {_display_path(output_path)}")
    if not os.access(parent, os.W_OK):
        raise GCDepthError(f"{label} directory is not writable: {_display_path(parent)}")
    return output_path


class AtomicTSVWriter:
    """Write a TSV to a temporary file and replace the destination only on success."""

    def __init__(self, destination: os.PathLike):
        self.destination = _prepare_output_path(destination, "Output-data")
        self._temporary_path: Optional[Path] = None
        self._handle: Optional[TextIO] = None
        self._committed = False

    def __enter__(self) -> "AtomicTSVWriter":
        temporary = tempfile.NamedTemporaryFile(
            mode="wt",
            encoding="utf-8",
            newline="",
            prefix=f".{self.destination.name}.",
            suffix=".tmp",
            dir=str(self.destination.parent),
            delete=False,
        )
        self._temporary_path = Path(temporary.name)
        self._handle = temporary
        self._handle.write("gc\tdepth\n")
        return self

    def write(self, gc_value: float, depth_value: float) -> None:
        if self._handle is None:
            raise RuntimeError("AtomicTSVWriter is not open")
        self._handle.write(f"{gc_value:.4f}\t{depth_value:.4f}\n")

    def commit(self) -> None:
        if self._handle is None or self._temporary_path is None:
            raise RuntimeError("AtomicTSVWriter is not open")
        self._handle.flush()
        os.fsync(self._handle.fileno())
        self._handle.close()
        self._handle = None
        os.replace(self._temporary_path, self.destination)
        self._committed = True

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None
        if not self._committed and self._temporary_path is not None:
            try:
                self._temporary_path.unlink()
            except FileNotFoundError:
                pass


class DiskPairStore:
    """Append float32 GC-depth pairs to disk and expose them as a NumPy memmap."""

    def __init__(self, directory: os.PathLike):
        self.path = Path(directory) / "plot-pairs.float32"
        self._handle = self.path.open("wb")
        self._buffer = array("f")
        self.count = 0
        self._closed = False

    def append(self, gc_value: float, depth_value: float) -> None:
        self._buffer.append(gc_value)
        self._buffer.append(depth_value)
        self.count += 1
        if len(self._buffer) >= PAIR_BUFFER_FLOATS:
            self._flush_buffer()

    def append_many(self, gc_values: np.ndarray, depth_values: np.ndarray) -> None:
        if len(gc_values) != len(depth_values):
            raise RuntimeError("GC and depth batch lengths differ")
        if len(gc_values) == 0:
            return
        self._flush_buffer()
        self._handle.flush()
        pairs = np.empty((len(gc_values), 2), dtype=np.float32)
        pairs[:, 0] = gc_values
        pairs[:, 1] = depth_values
        pairs.tofile(self._handle)
        self.count += len(gc_values)

    def _flush_buffer(self) -> None:
        if self._buffer:
            self._buffer.tofile(self._handle)
            self._buffer = array("f")

    def close(self) -> None:
        if not self._closed:
            self._flush_buffer()
            self._handle.flush()
            self._handle.close()
            self._closed = True

    def memmap(self) -> np.memmap:
        self.close()
        if self.count == 0:
            raise GCDepthError("No valid positive-depth data points are available for plotting.")
        return np.memmap(self.path, dtype=np.float32, mode="r", shape=(self.count, 2))

    def __enter__(self) -> "DiskPairStore":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


class DepthValueWriter:
    """Write compact float32 depth values to a temporary binary file in batches."""

    def __init__(self, path: os.PathLike):
        self.path = Path(path)
        self._handle = self.path.open("wb")
        self._depths = array("f")
        self.count = 0
        self._closed = False

    def append(self, depth: float) -> None:
        self._depths.append(depth)
        self.count += 1
        if len(self._depths) >= DEPTH_RECORD_BUFFER:
            self._flush()

    def _flush(self) -> None:
        if not self._depths:
            return
        self._depths.tofile(self._handle)
        self._depths = array("f")

    def close(self) -> None:
        if not self._closed:
            self._flush()
            self._handle.flush()
            self._handle.close()
            self._closed = True

    def __enter__(self) -> "DepthValueWriter":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


class DepthIndex:
    """A compact per-sequence SQLite index over disk-mapped PanDepth records."""

    def __init__(
        self,
        connection: sqlite3.Connection,
        records_path: os.PathLike,
        window_count: int,
        sequence_count: int,
    ):
        self.connection = connection
        self.records_path = Path(records_path)
        self.window_count = window_count
        self.sequence_count = sequence_count
        self._records: Optional[np.memmap] = None

    def _record_map(self) -> np.memmap:
        if self._records is None:
            self._records = np.memmap(
                self.records_path,
                dtype=np.float32,
                mode="r",
                shape=(self.window_count,),
            )
        return self._records

    def data_for_sequence(
        self, sequence_name: str
    ) -> Optional[Tuple[np.ndarray, int]]:
        row = self.connection.execute(
            "SELECT record_offset, record_count, sequence_length "
            "FROM sequences WHERE name = ?",
            (sequence_name,),
        ).fetchone()
        if row is None:
            return None

        record_offset = int(row[0])
        record_count = int(row[1])
        sequence_length = int(row[2])
        values = self._record_map()[record_offset : record_offset + record_count]
        return values, sequence_length

    def close(self) -> None:
        if self._records is not None:
            mmap_handle = getattr(self._records, "_mmap", None)
            self._records = None
            if mmap_handle is not None:
                mmap_handle.close()
        self.connection.close()


def _configure_sqlite(connection: sqlite3.Connection) -> None:
    """Configure the small temporary sequence index for low overhead."""
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-8192")
    connection.execute("PRAGMA locking_mode=EXCLUSIVE")


def _create_depth_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE TABLE sequences (
            name TEXT PRIMARY KEY,
            record_offset INTEGER NOT NULL,
            record_count INTEGER NOT NULL,
            sequence_length INTEGER NOT NULL
        ) WITHOUT ROWID;
        CREATE TABLE fasta_seen (
            name TEXT PRIMARY KEY
        ) WITHOUT ROWID;
        """
    )


def _parse_pandepth_row(
    line: str,
    line_number: int,
    window_size: int,
) -> Tuple[str, int, int, float]:
    """Parse and validate one PanDepth data row."""
    parts = line.rstrip("\r\n").split("\t")
    if len(parts) < 8:
        raise GCDepthError(
            f"Malformed PanDepth row at line {line_number}: expected at least 8 tab-separated columns, "
            f"found {len(parts)}."
        )

    sequence_name = parts[0].strip()
    if not sequence_name:
        raise GCDepthError(f"Empty sequence name in PanDepth file at line {line_number}.")

    try:
        start = int(parts[1])
        end = int(parts[2])
        length = int(parts[3])
        mean_depth = float(parts[7])
    except ValueError as exc:
        raise GCDepthError(
            f"Invalid numeric value in PanDepth file at line {line_number}."
        ) from exc

    if start < 1:
        raise GCDepthError(f"PanDepth start must be at least 1 at line {line_number}.")
    if end < start:
        raise GCDepthError(f"PanDepth end is smaller than start at line {line_number}.")
    if length <= 0:
        raise GCDepthError(f"PanDepth window length must be positive at line {line_number}.")
    if end - start + 1 != length:
        raise GCDepthError(
            f"Inconsistent PanDepth coordinates at line {line_number}: "
            f"end - start + 1 is not equal to Length."
        )
    if length > window_size:
        raise GCDepthError(
            f"PanDepth window length {length} exceeds --window {window_size} at line {line_number}. "
            "Use the same window size that was supplied to PanDepth."
        )
    if (start - 1) % window_size != 0:
        raise GCDepthError(
            f"PanDepth start {start} is not aligned to --window {window_size} at line {line_number}. "
            "Use the same window size that was supplied to PanDepth."
        )
    if not math.isfinite(mean_depth) or mean_depth < 0:
        raise GCDepthError(
            f"PanDepth MeanDepth must be a finite non-negative number at line {line_number}."
        )
    if mean_depth > FLOAT32_MAX:
        raise GCDepthError(
            f"PanDepth MeanDepth exceeds the supported float32 range at line {line_number}."
        )

    return sequence_name, start, length, mean_depth


def build_depth_index(
    pandepth_path: os.PathLike,
    database_path: os.PathLike,
    records_path: os.PathLike,
    window_size: int,
) -> DepthIndex:
    """Stream grouped PanDepth rows into compact binary storage plus a sequence index."""
    try:
        connection = sqlite3.connect(str(database_path))
    except sqlite3.Error as exc:
        raise GCDepthError(f"Cannot create temporary depth index: {exc}") from exc

    try:
        _configure_sqlite(connection)
        _create_depth_schema(connection)
        valid_rows = 0
        current_sequence: Optional[str] = None
        current_offset = 0
        current_count = 0
        previous_start = 0
        previous_length = 0

        def finish_sequence() -> None:
            nonlocal current_sequence, current_offset, current_count
            if current_sequence is None:
                return
            try:
                sequence_length = previous_start + previous_length - 1
                connection.execute(
                    "INSERT INTO sequences("
                    "name, record_offset, record_count, sequence_length"
                    ") VALUES (?, ?, ?, ?)",
                    (current_sequence, current_offset, current_count, sequence_length),
                )
            except sqlite3.IntegrityError as exc:
                raise GCDepthError(
                    f"PanDepth sequence '{current_sequence}' occurs in more than one block. "
                    "PanDepth windows must be grouped by sequence name."
                ) from exc

        with DepthValueWriter(records_path) as record_writer:
            with open_file(pandepth_path) as handle:
                for line_number, line in enumerate(handle, start=1):
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue

                    sequence_name, start, length, mean_depth = _parse_pandepth_row(
                        line, line_number, window_size
                    )

                    if sequence_name != current_sequence:
                        finish_sequence()
                        current_sequence = sequence_name
                        current_offset = record_writer.count
                        current_count = 0
                        previous_start = 0
                        previous_length = 0
                        if start != 1:
                            raise GCDepthError(
                                f"The first PanDepth window for '{sequence_name}' starts at {start}, "
                                f"but GC-depth expected 1 at line {line_number}."
                            )
                    else:
                        if previous_length != window_size:
                            raise GCDepthError(
                                f"A short PanDepth window for '{sequence_name}' appears before the final "
                                f"record at line {line_number}."
                            )
                        expected_start = previous_start + window_size
                        if start == previous_start:
                            raise GCDepthError(
                                f"Duplicate PanDepth window for {sequence_name}:{start} at line {line_number}."
                            )
                        if start != expected_start:
                            raise GCDepthError(
                                f"Non-contiguous PanDepth windows for '{sequence_name}' at line {line_number}: "
                                f"expected start {expected_start}, found {start}."
                            )

                    record_writer.append(mean_depth)
                    current_count += 1
                    valid_rows += 1
                    previous_start = start
                    previous_length = length

            finish_sequence()

        if valid_rows == 0:
            raise GCDepthError(
                "The PanDepth file contains no data rows. Expected tab-separated window records."
            )

        connection.commit()
        sequence_count = int(
            connection.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        )
        return DepthIndex(
            connection, records_path, valid_rows, sequence_count
        )

    except (UnicodeError, gzip.BadGzipFile, EOFError) as exc:
        connection.close()
        raise GCDepthError(f"Cannot decode or decompress PanDepth file: {exc}") from exc
    except sqlite3.Error as exc:
        connection.close()
        raise GCDepthError(f"Temporary depth index error: {exc}") from exc
    except Exception:
        connection.close()
        raise


def _register_fasta_name(connection: sqlite3.Connection, sequence_name: str) -> None:
    """Register a FASTA ID and fail clearly on duplicate IDs."""
    try:
        connection.execute("INSERT INTO fasta_seen(name) VALUES (?)", (sequence_name,))
    except sqlite3.IntegrityError as exc:
        raise GCDepthError(
            f"Duplicate FASTA sequence identifier: {sequence_name}. "
            "GC-depth uses the first word after '>' as the identifier."
        ) from exc


def compute_gc_windows_streaming(
    fasta_path: os.PathLike,
    depth_index: DepthIndex,
    window_size: int,
    pair_store: DiskPairStore,
    output_writer: Optional[AtomicTSVWriter] = None,
) -> Dict[str, object]:
    """Stream FASTA windows, merge against disk-backed depth rows, and store plot pairs on disk."""
    fasta_sequences = 0
    fasta_windows = 0
    matched_windows = 0
    zero_depth_windows = 0
    unmatched_fasta_windows = 0
    unmatched_fasta_sequences = 0
    unmatched_fasta_examples = []

    current_name: Optional[str] = None
    buffer = bytearray()
    sequence_batch = bytearray()
    depth_values: Optional[np.ndarray] = None
    depth_count = 0
    expected_sequence_length: Optional[int] = None
    current_sequence_windows = 0
    current_sequence_bases = 0
    saw_nonblank_content = False

    def start_sequence(header_line: str, line_number: int) -> None:
        nonlocal current_name, buffer, sequence_batch, depth_values
        nonlocal depth_count, expected_sequence_length
        nonlocal current_sequence_windows, current_sequence_bases
        nonlocal fasta_sequences, unmatched_fasta_sequences

        header = header_line[1:].strip()
        if not header:
            raise GCDepthError(f"Empty FASTA header at line {line_number}.")
        sequence_name = header.split()[0]
        _register_fasta_name(depth_index.connection, sequence_name)

        current_name = sequence_name
        buffer = bytearray()
        sequence_batch = bytearray()
        depth_data = depth_index.data_for_sequence(sequence_name)
        if depth_data is None:
            depth_values = None
            depth_count = 0
            expected_sequence_length = None
        else:
            depth_values, expected_sequence_length = depth_data
            depth_count = len(depth_values)
        current_sequence_windows = 0
        current_sequence_bases = 0
        fasta_sequences += 1

        if depth_count == 0:
            unmatched_fasta_sequences += 1
            if len(unmatched_fasta_examples) < 5:
                unmatched_fasta_examples.append(sequence_name)

    def record_gc_batch(gc_values: np.ndarray) -> None:
        nonlocal fasta_windows, matched_windows
        nonlocal zero_depth_windows, unmatched_fasta_windows, current_sequence_windows

        batch_count = len(gc_values)
        if batch_count == 0:
            return
        first_window_index = current_sequence_windows
        fasta_windows += batch_count
        current_sequence_windows += batch_count

        if depth_values is None or first_window_index >= depth_count:
            unmatched_fasta_windows += batch_count
            return

        matched_count = min(batch_count, depth_count - first_window_index)
        matched_gc = np.asarray(gc_values[:matched_count], dtype=np.float32)
        matched_depth = np.asarray(
            depth_values[first_window_index : first_window_index + matched_count],
            dtype=np.float32,
        )
        matched_windows += matched_count
        unmatched_fasta_windows += batch_count - matched_count

        if output_writer is not None:
            for gc_value, depth_value in zip(matched_gc, matched_depth):
                output_writer.write(float(gc_value), float(depth_value))

        positive_mask = matched_depth > 0
        positive_count = int(np.count_nonzero(positive_mask))
        zero_depth_windows += matched_count - positive_count
        if positive_count:
            pair_store.append_many(
                matched_gc[positive_mask], matched_depth[positive_mask]
            )

    def process_sequence_batch() -> None:
        nonlocal buffer, sequence_batch
        if not sequence_batch:
            return

        if buffer:
            combined = bytes(buffer) + bytes(sequence_batch)
        else:
            combined = bytes(sequence_batch)
        sequence_batch.clear()

        complete_length = (len(combined) // window_size) * window_size
        if complete_length:
            base_values = np.frombuffer(
                combined, dtype=np.uint8, count=complete_length
            ).reshape(-1, window_size)
            gc_counts = GC_LOOKUP[base_values].sum(axis=1, dtype=np.uint64)
            gc_values = (
                gc_counts.astype(np.float64) / float(window_size) * 100.0
            ).astype(np.float32)
            record_gc_batch(gc_values)

        buffer = bytearray(combined[complete_length:])

    def finish_sequence() -> None:
        nonlocal buffer
        if current_name is None:
            return
        if current_sequence_bases == 0:
            raise GCDepthError(f"FASTA sequence has no bases: {current_name}")
        if (
            expected_sequence_length is not None
            and current_sequence_bases != expected_sequence_length
        ):
            raise GCDepthError(
                f"Sequence length mismatch for '{current_name}': FASTA has "
                f"{current_sequence_bases:,} bases, but PanDepth describes "
                f"{expected_sequence_length:,} bases."
            )
        process_sequence_batch()
        if buffer:
            gc_count = buffer.count(ord("G")) + buffer.count(ord("C"))
            gc_value = np.array(
                [gc_count / len(buffer) * 100.0], dtype=np.float32
            )
            record_gc_batch(gc_value)
            buffer.clear()

    def accept_header(header_text: str, line_number: int) -> None:
        finish_sequence()
        start_sequence(">" + header_text, line_number)

    def accept_sequence(sequence_bytes: bytes, line_number: int) -> None:
        nonlocal current_sequence_bases
        if current_name is None:
            raise GCDepthError(
                f"FASTA sequence data appears before the first header at line {line_number}."
            )
        if not sequence_bytes.isascii():
            raise GCDepthError(
                f"Non-ASCII character in FASTA sequence at line {line_number}."
            )
        if INVALID_SEQUENCE_BYTES.search(sequence_bytes):
            raise GCDepthError(
                f"Whitespace, control character, or misplaced '>' in FASTA sequence "
                f"at line {line_number}."
            )

        uppercase = sequence_bytes.upper()
        current_sequence_bases += len(uppercase)
        sequence_batch.extend(uppercase)
        if len(sequence_batch) >= 1_048_576:
            process_sequence_batch()

    saw_nonblank_content = _stream_fasta_events(
        fasta_path, accept_header, accept_sequence
    )
    finish_sequence()

    if not saw_nonblank_content:
        raise GCDepthError("The FASTA file is empty.")
    if fasta_sequences == 0:
        raise GCDepthError("The FASTA file contains no valid sequence headers.")

    depth_index.connection.commit()
    total_depth_windows = depth_index.window_count
    unmatched_depth_windows = total_depth_windows - matched_windows
    missing_depth_sequences = int(
        depth_index.connection.execute(
            "SELECT COUNT(*) FROM sequences AS s "
            "LEFT JOIN fasta_seen AS f ON f.name = s.name "
            "WHERE f.name IS NULL"
        ).fetchone()[0]
    )
    missing_depth_examples = [
        str(row[0])
        for row in depth_index.connection.execute(
            "SELECT s.name FROM sequences AS s "
            "LEFT JOIN fasta_seen AS f ON f.name = s.name "
            "WHERE f.name IS NULL ORDER BY s.name LIMIT 5"
        )
    ]

    if matched_windows == 0:
        details = []
        if unmatched_fasta_examples:
            details.append("FASTA-only examples: " + ", ".join(unmatched_fasta_examples))
        if missing_depth_examples:
            details.append("PanDepth-only examples: " + ", ".join(missing_depth_examples))
        suffix = f" {'; '.join(details)}." if details else ""
        raise GCDepthError(
            "No windows matched between the FASTA and PanDepth inputs. "
            "Check sequence identifiers and ensure --window matches the PanDepth window size."
            + suffix
        )

    return {
        "fasta_sequences": fasta_sequences,
        "fasta_windows": fasta_windows,
        "matched_windows": matched_windows,
        "zero_depth_windows": zero_depth_windows,
        "unmatched_fasta_windows": unmatched_fasta_windows,
        "unmatched_fasta_sequences": unmatched_fasta_sequences,
        "unmatched_fasta_examples": unmatched_fasta_examples,
        "unmatched_depth_windows": unmatched_depth_windows,
        "missing_depth_sequences": missing_depth_sequences,
        "missing_depth_examples": missing_depth_examples,
    }


def load_combined_streaming(path: os.PathLike, pair_store: DiskPairStore) -> Dict[str, int]:
    """Validate and stream a combined TSV into disk-backed plotting storage."""
    data_rows = 0
    zero_depth_rows = 0
    saw_header = False

    try:
        with open_file(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped:
                    continue

                parts = line.rstrip("\r\n").split("\t")
                if not saw_header:
                    saw_header = True
                    if len(parts) >= 2 and parts[0].strip().lower() == "gc" and parts[1].strip().lower() == "depth":
                        continue

                if len(parts) < 2:
                    raise GCDepthError(
                        f"Malformed combined TSV row at line {line_number}: expected two tab-separated columns."
                    )
                try:
                    gc_value = float(parts[0])
                    depth_value = float(parts[1])
                except ValueError as exc:
                    raise GCDepthError(
                        f"Invalid numeric value in combined TSV at line {line_number}."
                    ) from exc

                if not math.isfinite(gc_value) or not 0 <= gc_value <= 100:
                    raise GCDepthError(
                        f"GC value must be between 0 and 100 at line {line_number}."
                    )
                if not math.isfinite(depth_value) or depth_value < 0:
                    raise GCDepthError(
                        f"Depth must be a finite non-negative number at line {line_number}."
                    )
                if depth_value > FLOAT32_MAX:
                    raise GCDepthError(
                        f"Depth exceeds the supported float32 range at line {line_number}."
                    )

                data_rows += 1
                if depth_value > 0:
                    pair_store.append(gc_value, depth_value)
                else:
                    zero_depth_rows += 1

    except (UnicodeError, gzip.BadGzipFile, EOFError) as exc:
        raise GCDepthError(f"Cannot decode or decompress combined TSV: {exc}") from exc

    if data_rows == 0:
        raise GCDepthError("The combined TSV contains no data rows.")

    return {"data_rows": data_rows, "zero_depth_rows": zero_depth_rows}


def _iter_pair_chunks(pairs: np.memmap, chunk_size: int = PLOT_CHUNK_SIZE):
    """Yield moderate-size contiguous pair chunks from a disk-backed array."""
    total = pairs.shape[0]
    for start in range(0, total, chunk_size):
        yield pairs[start : min(start + chunk_size, total)]


def _stream_histogram2d(
    pairs: np.memmap,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> np.ndarray:
    histogram = np.zeros((len(x_edges) - 1, len(y_edges) - 1), dtype=np.float64)
    for chunk in _iter_pair_chunks(pairs):
        partial, _, _ = np.histogram2d(
            chunk[:, 0], chunk[:, 1], bins=[x_edges, y_edges]
        )
        histogram += partial
    return histogram


def _stream_histogram(
    pairs: np.memmap,
    column: int,
    edges: np.ndarray,
) -> np.ndarray:
    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    for chunk in _iter_pair_chunks(pairs):
        partial, _ = np.histogram(chunk[:, column], bins=edges)
        counts += partial
    return counts


def _median_from_memmap(values: np.ndarray) -> float:
    """Calculate an exact median with one bounded in-memory float32 copy."""
    return float(np.median(np.asarray(values, dtype=np.float32)))


def _atomic_save_figure(fig: plt.Figure, output_path: Path) -> None:
    """Save the figure atomically so failed writes do not leave a partial plot."""
    suffix = output_path.suffix.lower()
    temporary = tempfile.NamedTemporaryFile(
        prefix=f".{output_path.stem}.",
        suffix=suffix,
        dir=str(output_path.parent),
        delete=False,
    )
    temporary_path = Path(temporary.name)
    temporary.close()
    try:
        fig.savefig(
            temporary_path,
            dpi=300,
            facecolor="white",
            edgecolor="none",
            bbox_inches="tight",
            pad_inches=0.1,
        )
        os.replace(temporary_path, output_path)
    except Exception:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass
        raise


def create_visualization(
    pairs: np.memmap,
    output_file: os.PathLike,
    log_depth: bool = False,
) -> Tuple[float, float, int]:
    """Generate the original GC-depth scatter plot using chunked disk-backed data access."""
    print("\n[Step 2] Creating visualization...")

    gc = pairs[:, 0]
    depth = pairs[:, 1]
    number_of_windows = int(pairs.shape[0])

    avg_gc = float(np.mean(gc, dtype=np.float64))
    avg_depth = float(np.mean(depth, dtype=np.float64))
    median_gc = _median_from_memmap(gc)
    median_depth = _median_from_memmap(depth)

    print(f"  GC average: {avg_gc:.2f}%, median: {median_gc:.2f}%")
    print(f"  Depth average: {avg_depth:.2f}X, median: {median_depth:.2f}X")
    print(f"  Total windows: {number_of_windows:,}")

    plt.style.use("seaborn-v0_8-whitegrid")
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
            "font.size": 14,
            "axes.labelsize": 16,
            "axes.titlesize": 14,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
            "axes.linewidth": 0.8,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "grid.alpha": 0.5,
            "grid.linestyle": "--",
            "grid.linewidth": 0.5,
            "grid.color": "#b0b0b0",
        }
    )

    colors = {
        "histogram": "#b8b8b8",
        "histogram_edge": "#909090",
        "average_line": "#DC143C",
        "grid": "#c0c0c0",
    }
    cmap = LinearSegmentedColormap.from_list(
        "gc_depth", ["#440154", "#21918c", "#fde725"]
    )

    gc_min = float(np.min(gc))
    gc_max = float(np.max(gc))
    depth_min_value = float(np.min(depth))
    gc_range_min = gc_min - 5
    gc_range_max = gc_max + 5
    depth_max = float(np.quantile(np.asarray(depth, dtype=np.float32), 0.99)) * 1.1
    if not math.isfinite(depth_max) or depth_max <= 0:
        raise GCDepthError("Could not determine a valid depth display range.")

    fig = plt.figure(figsize=(11, 10), facecolor="white")
    gs = GridSpec(
        4,
        4,
        figure=fig,
        height_ratios=[1.2, 1.2, 1.2, 1.2],
        width_ratios=[1, 1, 1, 1.3],
        hspace=0.45,
        wspace=0.55,
    )

    ax_main = fig.add_subplot(gs[1:4, 0:3])

    print("  Calculating point density in chunks...")
    x_edges = np.linspace(gc_range_min, gc_range_max, 100)
    if log_depth:
        density_depth_min = depth_min_value * 0.8
        y_edges = np.geomspace(density_depth_min, depth_max, 100)
    else:
        density_depth_min = 0.0
        y_edges = np.linspace(0, depth_max, 100)

    hist2d = _stream_histogram2d(pairs, x_edges, y_edges)
    hist2d_smooth = gaussian_filter(hist2d.T, sigma=1.5)

    density_min = math.inf
    density_max = -math.inf
    for chunk in _iter_pair_chunks(pairs):
        x_idx = np.clip(np.digitize(chunk[:, 0], x_edges) - 1, 0, len(x_edges) - 2)
        y_idx = np.clip(np.digitize(chunk[:, 1], y_edges) - 1, 0, len(y_edges) - 2)
        chunk_density = hist2d_smooth[y_idx, x_idx]
        density_min = min(density_min, float(np.min(chunk_density)))
        density_max = max(density_max, float(np.max(chunk_density)))

    density_span = density_max - density_min
    for chunk in _iter_pair_chunks(pairs):
        x_idx = np.clip(np.digitize(chunk[:, 0], x_edges) - 1, 0, len(x_edges) - 2)
        y_idx = np.clip(np.digitize(chunk[:, 1], y_edges) - 1, 0, len(y_edges) - 2)
        chunk_density = hist2d_smooth[y_idx, x_idx]
        if density_span > 0:
            density_norm = (chunk_density - density_min) / density_span
        else:
            density_norm = np.zeros_like(chunk_density)

        ax_main.scatter(
            chunk[:, 0],
            chunk[:, 1],
            c=density_norm,
            s=15,
            alpha=0.7,
            cmap=cmap,
            edgecolors="none",
            rasterized=True,
            vmin=0,
            vmax=1,
        )

    ax_main.axhline(
        y=median_depth,
        color=colors["average_line"],
        linestyle="--",
        linewidth=1.5,
        alpha=0.8,
    )
    ax_main.axvline(
        x=median_gc,
        color=colors["average_line"],
        linestyle="--",
        linewidth=1.5,
        alpha=0.8,
    )
    ax_main.set_xlabel(
        f"GC % (Average : {avg_gc:.2f} %)", fontsize=13, fontweight="bold"
    )
    ax_main.set_ylabel(
        f"Depth (Average : {avg_depth:.2f} X)", fontsize=13, fontweight="bold"
    )
    ax_main.set_xlim(gc_range_min, gc_range_max)
    if log_depth:
        ax_main.set_yscale("log")
        ax_main.set_ylim(density_depth_min, depth_max)
    else:
        ax_main.set_ylim(0, depth_max)
    ax_main.grid(
        True,
        linestyle="--",
        alpha=0.5,
        color=colors["grid"],
        linewidth=0.5,
    )
    ax_main.set_axisbelow(True)
    for spine in ax_main.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax_main.tick_params(
        axis="both",
        which="both",
        direction="out",
        top=False,
        bottom=True,
        left=True,
        right=False,
        length=5,
        width=0.8,
    )

    ax_gc = fig.add_subplot(gs[0, 0:3])
    gc_hist_edges = np.linspace(gc_range_min, gc_range_max, 61)
    gc_counts = _stream_histogram(pairs, 0, gc_hist_edges)
    ax_gc.bar(
        gc_hist_edges[:-1],
        gc_counts,
        width=np.diff(gc_hist_edges),
        align="edge",
        color=colors["histogram"],
        edgecolor=colors["histogram_edge"],
        linewidth=0.3,
        alpha=0.9,
    )
    ax_gc.axvline(
        x=median_gc,
        color=colors["average_line"],
        linestyle="--",
        linewidth=1.5,
        alpha=0.8,
    )
    ax_gc.set_ylabel("Numbers", fontsize=13, fontweight="bold")
    ax_gc.set_xlabel("")
    ax_gc.set_xlim(gc_range_min, gc_range_max)
    ax_gc.xaxis.set_tick_params(labelbottom=True)
    ax_gc.grid(
        True,
        linestyle="--",
        alpha=0.5,
        color=colors["grid"],
        linewidth=0.5,
    )
    ax_gc.set_axisbelow(True)
    for spine in ax_gc.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax_gc.tick_params(
        axis="both",
        which="both",
        direction="out",
        top=False,
        bottom=True,
        left=True,
        right=False,
        length=5,
        width=0.8,
    )

    ax_depth = fig.add_subplot(gs[1:4, 3])
    depth_range = (
        (density_depth_min, depth_max) if log_depth else (0.0, depth_max)
    )
    depth_hist_edges = (
        np.geomspace(depth_range[0], depth_range[1], 61)
        if log_depth
        else np.linspace(depth_range[0], depth_range[1], 61)
    )
    depth_counts = _stream_histogram(pairs, 1, depth_hist_edges)
    ax_depth.barh(
        depth_hist_edges[:-1],
        depth_counts,
        height=np.diff(depth_hist_edges),
        align="edge",
        color=colors["histogram"],
        edgecolor=colors["histogram_edge"],
        linewidth=0.3,
        alpha=0.9,
    )
    ax_depth.axhline(
        y=median_depth,
        color=colors["average_line"],
        linestyle="--",
        linewidth=1.5,
        alpha=0.8,
    )
    ax_depth.set_xlabel("Numbers", fontsize=13, fontweight="bold")
    ax_depth.tick_params(axis="x", labelsize=11)
    ax_depth.set_ylabel("")
    if log_depth:
        ax_depth.set_yscale("log")
    ax_depth.set_ylim(depth_range)
    ax_depth.yaxis.set_tick_params(labelleft=True)
    ax_depth.grid(
        True,
        linestyle="--",
        alpha=0.5,
        color=colors["grid"],
        linewidth=0.5,
    )
    ax_depth.set_axisbelow(True)
    for spine in ax_depth.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax_depth.tick_params(
        axis="both",
        which="both",
        direction="out",
        top=False,
        bottom=True,
        left=True,
        right=False,
        length=5,
        width=0.8,
    )

    fig.subplots_adjust(
        left=0.10,
        right=0.95,
        bottom=0.08,
        top=0.95,
        hspace=0.5,
        wspace=0.5,
    )

    output_path = _prepare_output_path(output_file, "Plot output")
    print(f"  Saving plot to: {output_path}")
    try:
        _atomic_save_figure(fig, output_path)
    except OSError as exc:
        raise GCDepthError(f"Could not save plot to {output_path}: {exc}") from exc
    finally:
        plt.close(fig)

    return median_gc, median_depth, number_of_windows


def _print_matching_warnings(stats: Dict[str, object]) -> None:
    unmatched_fasta_windows = int(stats["unmatched_fasta_windows"])
    unmatched_depth_windows = int(stats["unmatched_depth_windows"])
    zero_depth_windows = int(stats["zero_depth_windows"])

    if unmatched_fasta_windows or unmatched_depth_windows:
        print("\n  Warning: not every input window was matched.")
        print(f"    FASTA windows without PanDepth data: {unmatched_fasta_windows:,}")
        print(f"    PanDepth windows without FASTA data: {unmatched_depth_windows:,}")

        fasta_examples = stats.get("unmatched_fasta_examples") or []
        depth_examples = stats.get("missing_depth_examples") or []
        if fasta_examples:
            print("    FASTA-only sequence examples: " + ", ".join(fasta_examples))
        if depth_examples:
            print("    PanDepth-only sequence examples: " + ", ".join(depth_examples))

    if zero_depth_windows:
        print(
            f"  Zero-depth windows excluded from the plot: {zero_depth_windows:,} "
            "(retained in --output-data when requested)"
        )


def main() -> None:
    parser = DocstringHelpParser(
        prog="gc_depth",
        description="Compute and visualize GC content versus sequencing depth per genomic window.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "fasta",
        nargs="?",
        help="Genome FASTA file, chromosome or contig level (.fa/.fasta, plain or gzipped)",
    )
    parser.add_argument(
        "pandepth",
        nargs="?",
        help="Pandepth windowed depth file (.win.stat.gz or plain text with the same format)",
    )
    parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=1000,
        help="Window size in bp, must match pandepth -w value (default: 1000)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="gc-depth.png",
        help="Output plot file (.png or .pdf, default: gc-depth.png)",
    )
    parser.add_argument(
        "--log-depth", action="store_true", help="Use logarithmic scale for the depth axis"
    )
    parser.add_argument(
        "--plot-only",
        metavar="TSV",
        help="Skip processing, re-plot from an existing combined TSV (from --output-data)",
    )
    parser.add_argument(
        "--output-data",
        metavar="FILE",
        help="Save merged GC and depth data to this TSV file (can be reused with --plot-only)",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    args = parser.parse_args()

    if args.window <= 0:
        parser.error("--window must be a positive integer")

    output_path = Path(args.output).expanduser()
    if output_path.suffix.lower() not in (".png", ".pdf"):
        parser.error("--output must use a .png or .pdf extension")

    if args.plot_only and (args.fasta or args.pandepth):
        parser.error("Do not provide fasta or pandepth positional arguments with --plot-only")
    if args.plot_only and args.output_data:
        parser.error("--output-data cannot be used together with --plot-only")
    if not args.plot_only and (not args.fasta or not args.pandepth):
        parser.error("fasta and pandepth arguments are required unless using --plot-only")

    print("gc_depth")

    try:
        _prepare_output_path(output_path, "Plot output")
        if args.output_data:
            _prepare_output_path(args.output_data, "Output-data")

        input_paths = []
        if args.plot_only:
            input_paths.append(Path(args.plot_only).expanduser())
        else:
            if args.fasta:
                input_paths.append(Path(args.fasta).expanduser())
            if args.pandepth:
                input_paths.append(Path(args.pandepth).expanduser())
        resolved_inputs = {path.resolve(strict=False) for path in input_paths}
        if output_path.resolve(strict=False) in resolved_inputs:
            raise GCDepthError("Plot output must not overwrite an input file.")
        if args.output_data:
            data_path = Path(args.output_data).expanduser()
            if data_path.resolve(strict=False) in resolved_inputs:
                raise GCDepthError("--output-data must not overwrite an input file.")
            if data_path.resolve(strict=False) == output_path.resolve(strict=False):
                raise GCDepthError("--output-data and --output must use different paths.")

        with tempfile.TemporaryDirectory(prefix="gc_depth_") as temporary_directory:
            with DiskPairStore(temporary_directory) as pair_store:

                if args.plot_only:
                    print(f"\nPlot-only mode: {args.plot_only}")
                    plot_stats = load_combined_streaming(args.plot_only, pair_store)
                    if plot_stats["zero_depth_rows"]:
                        print(
                            f"  Zero-depth rows excluded from the plot: "
                            f"{plot_stats['zero_depth_rows']:,}"
                        )
                else:
                    print("\nParameters:")
                    print(f"  FASTA:       {args.fasta}")
                    print(f"  Depth file:  {args.pandepth}")
                    print(f"  Window size: {args.window} bp")
                    print(f"  Output:      {args.output}")
    
                    print("\n[Step 1] Processing data...")
                    print("  Streaming PanDepth data into temporary disk storage...")
                    database_path = Path(temporary_directory) / "depth-index.sqlite3"
                    records_path = Path(temporary_directory) / "depth-records.bin"
                    depth_index = build_depth_index(
                        args.pandepth, database_path, records_path, args.window
                    )
                    print(f"  Valid PanDepth windows: {depth_index.window_count:,}")
                    print(f"  PanDepth sequences: {depth_index.sequence_count:,}")
    
                    try:
                        print("  Streaming FASTA windows and matching depth records...")
                        if args.output_data:
                            with AtomicTSVWriter(args.output_data) as output_writer:
                                matching_stats = compute_gc_windows_streaming(
                                    args.fasta,
                                    depth_index,
                                    args.window,
                                    pair_store,
                                    output_writer,
                                )
                                output_writer.commit()
                            print(f"  Combined data saved to: {args.output_data}")
                        else:
                            matching_stats = compute_gc_windows_streaming(
                                args.fasta,
                                depth_index,
                                args.window,
                                pair_store,
                            )
                    finally:
                        depth_index.close()
    
                    print(
                        f"  Matched windows with GC and depth data: "
                        f"{matching_stats['matched_windows']:,}"
                    )
                    _print_matching_warnings(matching_stats)
    
                pairs = pair_store.memmap()
                median_gc, median_depth, number_of_windows = create_visualization(
                    pairs, output_path, args.log_depth
                )
                mmap_handle = getattr(pairs, "_mmap", None)
                del pairs
                if mmap_handle is not None:
                    mmap_handle.close()
        print("\nResults:")
        print(f"- Visualization output : {args.output}")
        print(f"- Median GC            : {median_gc:.2f}%")
        print(f"- Median Depth         : {median_depth:.2f}X")
        print(f"- Plotted windows      : {number_of_windows:,}")
        print("\nDone!")

    except GCDepthError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    except KeyboardInterrupt:
        print("\nError: Interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except MemoryError as exc:
        print(
            "Error: Not enough memory to finish plotting. Try a larger system or reduce the number of windows.",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc
    except (OSError, sqlite3.Error) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
