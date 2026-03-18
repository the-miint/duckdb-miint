#!/usr/bin/env python3
"""Generate mzXML test data files for read_mzxml unit tests.

Run once, check results into data/mzxml/. Not a build step.

mzXML binary encoding: big-endian IEEE 754, interleaved m/z-intensity pairs.
"""

import base64
import struct
import zlib
import os

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "mzxml")

MZXML_HEADER = '<?xml version="1.0" encoding="ISO-8859-1"?>\n'


def encode_peaks_64(mz_values, int_values, compress=False):
    """Encode m/z and intensity as big-endian 64-bit interleaved pairs."""
    data = b""
    for mz, intensity in zip(mz_values, int_values):
        data += struct.pack(">d", mz)  # big-endian double
        data += struct.pack(">d", intensity)
    if compress:
        data = zlib.compress(data)
    return base64.b64encode(data).decode("ascii")


def encode_peaks_32(mz_values, int_values):
    """Encode m/z and intensity as big-endian 32-bit interleaved pairs."""
    data = b""
    for mz, intensity in zip(mz_values, int_values):
        data += struct.pack(">f", mz)  # big-endian float
        data += struct.pack(">f", intensity)
    return base64.b64encode(data).decode("ascii")


def write_file(name, content):
    path = os.path.join(OUTPUT_DIR, name)
    with open(path, "w") as f:
        f.write(content)
    print(f"  wrote {path}")


def generate_basic_3spectra():
    """Flat layout: MS1 (centroid, positive), MS2 (CID), MS2 (profile, negative, HCD)."""
    mz = [100.0, 200.0, 300.0]
    intensities = [1000.0, 5000.0, 2000.0]
    peaks1 = encode_peaks_64(mz, intensities)

    mz2 = [195.0, 200.0, 205.0]
    int2 = [500.0, 3000.0, 800.0]
    peaks2 = encode_peaks_64(mz2, int2)

    mz3 = [290.0, 300.0, 310.0]
    int3 = [400.0, 2500.0, 600.0]
    peaks3 = encode_peaks_64(mz3, int3)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="3"
        retentionTime="PT90S"
        centroided="1"
        polarity="+"
        basePeakMz="200.0"
        basePeakIntensity="5000.0"
        totIonCurrent="8000.0"
        lowMz="100.0"
        highMz="300.0"
        filterLine="FTMS + p NSI Full ms [100.00-500.00]"
        startMz="100.0"
        endMz="500.0">
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none"
          compressedLen="0">{peaks1}</peaks>
  </scan>
  <scan num="2"
        msLevel="2"
        peaksCount="3"
        retentionTime="PT95S"
        centroided="1"
        polarity="+">
   <precursorMz precursorIntensity="5000.0"
                precursorCharge="2"
                activationMethod="CID"
                windowWideness="3.0"
                precursorScanNum="1">200.0</precursorMz>
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none"
          compressedLen="0">{peaks2}</peaks>
  </scan>
  <scan num="3"
        msLevel="2"
        peaksCount="3"
        retentionTime="PT108S"
        centroided="0"
        polarity="-">
   <precursorMz precursorIntensity="2000.0"
                precursorCharge="3"
                activationMethod="HCD"
                windowWideness="4.0"
                precursorScanNum="1">300.0</precursorMz>
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none"
          compressedLen="0">{peaks3}</peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("basic_3spectra.mzXML", content)


def generate_compressed():
    """MS1 with zlib compression, same values as basic MS1."""
    mz = [100.0, 200.0, 300.0]
    intensities = [1000.0, 5000.0, 2000.0]
    peaks = encode_peaks_64(mz, intensities, compress=True)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="3"
        retentionTime="PT90S"
        centroided="1"
        polarity="+"
        basePeakMz="200.0"
        basePeakIntensity="5000.0"
        totIonCurrent="8000.0"
        lowMz="100.0"
        highMz="300.0">
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="zlib">{peaks}</peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("compressed.mzXML", content)


def generate_float32():
    """MS1 with 32-bit precision."""
    mz = [100.0, 200.0, 300.0]
    intensities = [1000.0, 5000.0, 2000.0]
    peaks = encode_peaks_32(mz, intensities)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="3"
        retentionTime="PT90S"
        centroided="1"
        polarity="+"
        basePeakMz="200.0"
        basePeakIntensity="5000.0"
        totIonCurrent="8000.0"
        lowMz="100.0"
        highMz="300.0">
   <peaks precision="32"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none">{peaks}</peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("float32.mzXML", content)


def generate_missing_optional():
    """Minimal MS1: only msLevel, peaksCount, peaks. No compressionType attribute."""
    mz = [150.0]
    intensities = [999.0]
    peaks = encode_peaks_64(mz, intensities)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="1">
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int">{peaks}</peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("missing_optional.mzXML", content)


def generate_empty():
    """Valid msRun with zero scan elements."""
    content = (
        MZXML_HEADER
        + """\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
 </msRun>
</mzXML>
"""
    )
    write_file("empty.mzXML", content)


def generate_orphan_ms2():
    """Two MS2 scans at top level (flat), no MS1."""
    mz = [200.0, 300.0]
    intensities = [1000.0, 2000.0]
    peaks = encode_peaks_64(mz, intensities)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="2"
        peaksCount="2">
   <precursorMz precursorCharge="2"
                activationMethod="CID">200.0</precursorMz>
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none">{peaks}</peaks>
  </scan>
  <scan num="2"
        msLevel="2"
        peaksCount="2">
   <precursorMz precursorCharge="3"
                activationMethod="HCD">300.0</precursorMz>
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none">{peaks}</peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("orphan_ms2.mzXML", content)


def generate_nested_scans():
    """MS1 containing two nested MS2 scan children."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [1000.0, 5000.0, 2000.0]
    peaks1 = encode_peaks_64(mz1, int1)

    mz2 = [195.0, 200.0, 205.0]
    int2 = [500.0, 3000.0, 800.0]
    peaks2 = encode_peaks_64(mz2, int2)

    mz3 = [290.0, 300.0, 310.0]
    int3 = [400.0, 2500.0, 600.0]
    peaks3 = encode_peaks_64(mz3, int3)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="3"
        retentionTime="PT90S"
        centroided="1"
        polarity="+"
        basePeakMz="200.0"
        basePeakIntensity="5000.0"
        totIonCurrent="8000.0"
        lowMz="100.0"
        highMz="300.0">
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none">{peaks1}</peaks>
   <scan num="2"
         msLevel="2"
         peaksCount="3"
         retentionTime="PT95S">
    <precursorMz precursorIntensity="5000.0"
                 precursorCharge="2"
                 activationMethod="CID"
                 windowWideness="3.0">200.0</precursorMz>
    <peaks precision="64"
           byteOrder="network"
           contentType="m/z-int"
           compressionType="none">{peaks2}</peaks>
   </scan>
   <scan num="3"
         msLevel="2"
         peaksCount="3"
         retentionTime="PT100S">
    <precursorMz precursorIntensity="2000.0"
                 precursorCharge="3"
                 activationMethod="HCD"
                 windowWideness="4.0">300.0</precursorMz>
    <peaks precision="64"
           byteOrder="network"
           contentType="m/z-int"
           compressionType="none">{peaks3}</peaks>
   </scan>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("nested_scans.mzXML", content)


def generate_precursorscannum():
    """Flat layout: MS1 (scan 1), MS2 (scan 2, precursorScanNum=1)."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [1000.0, 5000.0, 2000.0]
    peaks1 = encode_peaks_64(mz1, int1)

    mz2 = [195.0, 200.0, 205.0]
    int2 = [500.0, 3000.0, 800.0]
    peaks2 = encode_peaks_64(mz2, int2)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="3"
        retentionTime="PT90S">
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none">{peaks1}</peaks>
  </scan>
  <scan num="2"
        msLevel="2"
        peaksCount="3"
        retentionTime="PT95S">
   <precursorMz precursorIntensity="5000.0"
                precursorCharge="2"
                activationMethod="CID"
                windowWideness="3.0"
                precursorScanNum="1">200.0</precursorMz>
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none">{peaks2}</peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("precursorscannum.mzXML", content)


def generate_zero_peaks():
    """MS1 with peaksCount=0 and empty peaks element."""
    content = (
        MZXML_HEADER
        + """\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="0"
        retentionTime="PT60S">
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none"></peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("zero_peaks.mzXML", content)


def generate_combined_rt():
    """MS1 with retentionTime=PT1M30S (combined minute+second = 1.5 min)."""
    mz = [150.0]
    intensities = [999.0]
    peaks = encode_peaks_64(mz, intensities)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="1"
        retentionTime="PT1M30S">
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none">{peaks}</peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("combined_rt.mzXML", content)


def generate_hours_rt():
    """MS1 with retentionTime=PT1H30M (hours + minutes = 90 min)."""
    mz = [150.0]
    intensities = [999.0]
    peaks = encode_peaks_64(mz, intensities)

    content = (
        MZXML_HEADER
        + f"""\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
 <msRun>
  <scan num="1"
        msLevel="1"
        peaksCount="1"
        retentionTime="PT1H30M">
   <peaks precision="64"
          byteOrder="network"
          contentType="m/z-int"
          compressionType="none">{peaks}</peaks>
  </scan>
 </msRun>
</mzXML>
"""
    )
    write_file("hours_rt.mzXML", content)


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("Generating mzXML test data...")
    generate_basic_3spectra()
    generate_compressed()
    generate_float32()
    generate_missing_optional()
    generate_empty()
    generate_orphan_ms2()
    generate_nested_scans()
    generate_precursorscannum()
    generate_zero_peaks()
    generate_combined_rt()
    generate_hours_rt()
    print("Done.")
