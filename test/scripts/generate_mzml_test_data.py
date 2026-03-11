#!/usr/bin/env python3
"""Generate synthetic mzML test data files for unit tests.

Each file tests a specific aspect of the mzML parser. The files are small
and self-contained, with known values for deterministic testing.
"""

import base64
import os
import struct
import zlib

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "mzml")


def encode_doubles_64(values):
    """Encode a list of doubles as base64, uncompressed, 64-bit LE."""
    raw = struct.pack(f"<{len(values)}d", *values)
    return base64.b64encode(raw).decode()


def encode_floats_32(values):
    """Encode a list of floats as base64, uncompressed, 32-bit LE."""
    raw = struct.pack(f"<{len(values)}f", *values)
    return base64.b64encode(raw).decode()


def encode_doubles_64_compressed(values):
    """Encode a list of doubles as base64, zlib-compressed, 64-bit LE."""
    raw = struct.pack(f"<{len(values)}d", *values)
    compressed = zlib.compress(raw)
    return base64.b64encode(compressed).decode()


def encode_floats_32_compressed(values):
    """Encode a list of floats as base64, zlib-compressed, 32-bit LE."""
    raw = struct.pack(f"<{len(values)}f", *values)
    compressed = zlib.compress(raw)
    return base64.b64encode(compressed).decode()


INDEXED_MZML_HEADER = """<?xml version="1.0" encoding="utf-8"?>
<indexedmzML xmlns="http://psi.hupo.org/ms/mzml"
             xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
             xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd">
<mzML xmlns="http://psi.hupo.org/ms/mzml">
<cvList count="2">
  <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>
  <cv id="UO" fullName="Unit Ontology" URI="http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/phenotype/unit.obo"/>
</cvList>
<fileDescription>
  <fileContent>
    <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum"/>
    <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum"/>
  </fileContent>
</fileDescription>
"""

INDEXED_MZML_FOOTER = """</mzML>
</indexedmzML>
"""

RUN_OPEN = '<run defaultInstrumentConfigurationRef="IC1">'
RUN_CLOSE = "</run>"


def make_binary_data_array(values, is_mz=True, precision=64, compressed=False):
    """Create a binaryDataArray XML element."""
    if compressed:
        if precision == 64:
            encoded = encode_doubles_64_compressed(values)
        else:
            encoded = encode_floats_32_compressed(values)
        compression_cv = '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression"/>'
    else:
        if precision == 64:
            encoded = encode_doubles_64(values)
        else:
            encoded = encode_floats_32(values)
        compression_cv = '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>'

    if precision == 64:
        precision_cv = '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
    else:
        precision_cv = '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"/>'

    if is_mz:
        array_cv = '<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>'
    else:
        array_cv = '<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>'

    return f"""        <binaryDataArray encodedLength="{len(encoded)}">
          {compression_cv}
          {precision_cv}
          {array_cv}
          <binary>{encoded}</binary>
        </binaryDataArray>"""


def make_spectrum(
    index,
    spectrum_id,
    ms_level,
    rt_value,
    rt_unit="minute",
    spectrum_type="centroid",
    polarity="positive",
    base_peak_mz=None,
    base_peak_intensity=None,
    tic=None,
    lowest_mz=None,
    highest_mz=None,
    mz_values=None,
    intensity_values=None,
    precursor_mz=None,
    precursor_charge=None,
    precursor_intensity=None,
    isolation_target=None,
    isolation_lower=None,
    isolation_upper=None,
    activation_method=None,
    collision_energy=None,
    filter_string=None,
    scan_window_lower=None,
    scan_window_upper=None,
    precision=64,
    compressed=False,
):
    """Build a spectrum XML element."""
    if mz_values is None:
        mz_values = []
    if intensity_values is None:
        intensity_values = []

    array_length = len(mz_values)

    lines = []
    lines.append(f'    <spectrum index="{index}" id="{spectrum_id}" defaultArrayLength="{array_length}">')
    lines.append(f'      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{ms_level}"/>')

    # spectrum type
    if spectrum_type == "centroid":
        lines.append('      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>')
    elif spectrum_type == "profile":
        lines.append('      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum"/>')

    # polarity
    if polarity == "positive":
        lines.append('      <cvParam cvRef="MS" accession="MS:1000130" name="positive scan"/>')
    elif polarity == "negative":
        lines.append('      <cvParam cvRef="MS" accession="MS:1000129" name="negative scan"/>')

    if base_peak_mz is not None:
        lines.append(f'      <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="{base_peak_mz}"/>')
    if base_peak_intensity is not None:
        lines.append(
            f'      <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="{base_peak_intensity}"/>'
        )
    if tic is not None:
        lines.append(f'      <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="{tic}"/>')
    if lowest_mz is not None:
        lines.append(
            f'      <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="{lowest_mz}"/>'
        )
    if highest_mz is not None:
        lines.append(
            f'      <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="{highest_mz}"/>'
        )
    if filter_string is not None:
        lines.append(f'      <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="{filter_string}"/>')

    # scanList
    lines.append('      <scanList count="1">')
    lines.append("        <scan>")
    if rt_unit == "minute":
        lines.append(
            f'          <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt_value}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>'
        )
    elif rt_unit == "second":
        lines.append(
            f'          <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt_value}" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>'
        )

    if scan_window_lower is not None or scan_window_upper is not None:
        lines.append('          <scanWindowList count="1">')
        lines.append("            <scanWindow>")
        if scan_window_lower is not None:
            lines.append(
                f'              <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="{scan_window_lower}"/>'
            )
        if scan_window_upper is not None:
            lines.append(
                f'              <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="{scan_window_upper}"/>'
            )
        lines.append("            </scanWindow>")
        lines.append("          </scanWindowList>")

    lines.append("        </scan>")
    lines.append("      </scanList>")

    # precursorList
    if precursor_mz is not None:
        lines.append('      <precursorList count="1">')
        lines.append("        <precursor>")
        if isolation_target is not None or isolation_lower is not None:
            lines.append("          <isolationWindow>")
            if isolation_target is not None:
                lines.append(
                    f'            <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{isolation_target}"/>'
                )
            if isolation_lower is not None:
                lines.append(
                    f'            <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{isolation_lower}"/>'
                )
            if isolation_upper is not None:
                lines.append(
                    f'            <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{isolation_upper}"/>'
                )
            lines.append("          </isolationWindow>")

        lines.append("          <selectedIonList>")
        lines.append("            <selectedIon>")
        lines.append(
            f'              <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{precursor_mz}"/>'
        )
        if precursor_charge is not None:
            lines.append(
                f'              <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{precursor_charge}"/>'
            )
        if precursor_intensity is not None:
            lines.append(
                f'              <cvParam cvRef="MS" accession="MS:1000042" name="peak intensity" value="{precursor_intensity}"/>'
            )
        lines.append("            </selectedIon>")
        lines.append("          </selectedIonList>")

        if activation_method is not None or collision_energy is not None:
            lines.append("          <activation>")
            if activation_method == "CID":
                lines.append(
                    '            <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation"/>'
                )
            elif activation_method == "HCD":
                lines.append(
                    '            <cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation"/>'
                )
            elif activation_method == "ETD":
                lines.append(
                    '            <cvParam cvRef="MS" accession="MS:1000598" name="electron transfer dissociation"/>'
                )
            if collision_energy is not None:
                lines.append(
                    f'            <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="{collision_energy}"/>'
                )
            lines.append("          </activation>")

        lines.append("        </precursor>")
        lines.append("      </precursorList>")

    # binaryDataArrayList
    lines.append(f'      <binaryDataArrayList count="2">')
    lines.append(make_binary_data_array(mz_values, is_mz=True, precision=precision, compressed=compressed))
    lines.append(make_binary_data_array(intensity_values, is_mz=False, precision=precision, compressed=compressed))
    lines.append("      </binaryDataArrayList>")
    lines.append("    </spectrum>")

    return "\n".join(lines)


def make_chromatogram(
    index,
    chrom_id,
    chrom_type,
    time_values,
    intensity_values,
    precursor_mz=None,
    product_mz=None,
    precision=64,
    compressed=False,
):
    """Build a chromatogram XML element."""
    array_length = len(time_values)
    lines = []
    lines.append(f'    <chromatogram index="{index}" id="{chrom_id}" defaultArrayLength="{array_length}">')

    # chromatogram type
    if chrom_type == "TIC":
        lines.append('      <cvParam cvRef="MS" accession="MS:1000235" name="total ion current chromatogram"/>')
    elif chrom_type == "BPC":
        lines.append('      <cvParam cvRef="MS" accession="MS:1000628" name="basepeak chromatogram"/>')
    elif chrom_type == "SRM":
        lines.append(
            '      <cvParam cvRef="MS" accession="MS:1001473" name="selected reaction monitoring chromatogram"/>'
        )
    elif chrom_type == "SIC":
        lines.append('      <cvParam cvRef="MS" accession="MS:1000627" name="selected ion current chromatogram"/>')

    # precursor
    if precursor_mz is not None:
        lines.append("      <precursor>")
        lines.append("        <isolationWindow>")
        lines.append(
            f'          <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{precursor_mz}"/>'
        )
        lines.append("        </isolationWindow>")
        lines.append("      </precursor>")

    # product
    if product_mz is not None:
        lines.append("      <product>")
        lines.append("        <isolationWindow>")
        lines.append(
            f'          <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{product_mz}"/>'
        )
        lines.append("        </isolationWindow>")
        lines.append("      </product>")

    # binary data arrays - time and intensity
    lines.append(f'      <binaryDataArrayList count="2">')

    # Time array
    if compressed:
        encoded_time = encode_doubles_64_compressed(time_values)
        compression_cv = '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression"/>'
    else:
        encoded_time = encode_doubles_64(time_values)
        compression_cv = '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>'

    if precision == 64:
        precision_cv = '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
    else:
        precision_cv = '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"/>'

    lines.append(f'        <binaryDataArray encodedLength="{len(encoded_time)}">')
    lines.append(f"          {compression_cv}")
    lines.append(f"          {precision_cv}")
    lines.append(
        '          <cvParam cvRef="MS" accession="MS:1000595" name="time array" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>'
    )
    lines.append(f"          <binary>{encoded_time}</binary>")
    lines.append("        </binaryDataArray>")

    # Intensity array
    lines.append(make_binary_data_array(intensity_values, is_mz=False, precision=precision, compressed=compressed))

    lines.append("      </binaryDataArrayList>")
    lines.append("    </chromatogram>")

    return "\n".join(lines)


def generate_basic_3spectra():
    """1 MS1 + 2 MS2 with various metadata. Uncompressed 64-bit."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [1000.0, 5000.0, 2000.0]

    mz2 = [150.0, 250.0]
    int2 = [3000.0, 1500.0]

    mz3 = [120.0, 220.0, 320.0, 420.0]
    int3 = [500.0, 2500.0, 1200.0, 800.0]

    spectra = []
    # Spectrum 0: MS1, centroid, positive, RT=1.5 min
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.5,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=200.0,
            base_peak_intensity=5000.0,
            tic=8000.0,
            lowest_mz=100.0,
            highest_mz=300.0,
            mz_values=mz1,
            intensity_values=int1,
            filter_string="FTMS + p NSI Full ms [100.00-500.00]",
            scan_window_lower=100.0,
            scan_window_upper=500.0,
        )
    )

    # Spectrum 1: MS2, centroid, positive, RT=1.6 min, CID, precursor 200.0
    spectra.append(
        make_spectrum(
            index=1,
            spectrum_id="scan=2",
            ms_level=2,
            rt_value=1.6,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=150.0,
            base_peak_intensity=3000.0,
            tic=4500.0,
            lowest_mz=150.0,
            highest_mz=250.0,
            mz_values=mz2,
            intensity_values=int2,
            precursor_mz=200.0,
            precursor_charge=2,
            precursor_intensity=5000.0,
            isolation_target=200.0,
            isolation_lower=1.5,
            isolation_upper=1.5,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    # Spectrum 2: MS2, profile, negative, RT=108 seconds (=1.8 min), HCD
    spectra.append(
        make_spectrum(
            index=2,
            spectrum_id="scan=3",
            ms_level=2,
            rt_value=108.0,
            rt_unit="second",
            spectrum_type="profile",
            polarity="negative",
            base_peak_mz=220.0,
            base_peak_intensity=2500.0,
            tic=5000.0,
            lowest_mz=120.0,
            highest_mz=420.0,
            mz_values=mz3,
            intensity_values=int3,
            precursor_mz=300.0,
            precursor_charge=3,
            precursor_intensity=2000.0,
            isolation_target=300.0,
            isolation_lower=2.0,
            isolation_upper=2.0,
            activation_method="HCD",
            collision_energy=30.0,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += f'  <spectrumList count="3" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "basic_3spectra.mzML"), "w") as f:
        f.write(content)


def generate_compressed():
    """Same 3 spectra but zlib-compressed binary arrays."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [1000.0, 5000.0, 2000.0]

    spectra = []
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.5,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=200.0,
            base_peak_intensity=5000.0,
            tic=8000.0,
            mz_values=mz1,
            intensity_values=int1,
            compressed=True,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += f'  <spectrumList count="1" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "compressed.mzML"), "w") as f:
        f.write(content)


def generate_float32():
    """32-bit float binary arrays."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [1000.0, 5000.0, 2000.0]

    spectra = []
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.5,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=200.0,
            base_peak_intensity=5000.0,
            tic=8000.0,
            mz_values=mz1,
            intensity_values=int1,
            precision=32,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += f'  <spectrumList count="1" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "float32.mzML"), "w") as f:
        f.write(content)


def generate_empty_spectra():
    """Valid mzML with spectrumList count=0."""
    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="0" defaultDataProcessingRef="dp">\n'
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "empty_spectra.mzML"), "w") as f:
        f.write(content)


def generate_missing_optional():
    """Only ms_level set, no polarity/base_peak/etc."""
    mz1 = [150.0]
    int1 = [1000.0]

    lines = []
    lines.append(f'    <spectrum index="0" id="scan=1" defaultArrayLength="1">')
    lines.append('      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>')
    lines.append('      <scanList count="1">')
    lines.append("        <scan>")
    lines.append("        </scan>")
    lines.append("      </scanList>")
    lines.append('      <binaryDataArrayList count="2">')
    lines.append(make_binary_data_array(mz1, is_mz=True))
    lines.append(make_binary_data_array(int1, is_mz=False))
    lines.append("      </binaryDataArrayList>")
    lines.append("    </spectrum>")

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="1" defaultDataProcessingRef="dp">\n'
    content += "\n".join(lines) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "missing_optional.mzML"), "w") as f:
        f.write(content)


def generate_with_chromatograms():
    """TIC + SRM + BPC chromatograms."""
    time1 = [0.5, 1.0, 1.5, 2.0]
    int1 = [10000.0, 25000.0, 30000.0, 15000.0]

    time2 = [0.5, 1.0, 1.5]
    int2 = [500.0, 1200.0, 800.0]

    time3 = [0.5, 1.0, 1.5, 2.0]
    int3 = [5000.0, 12000.0, 15000.0, 8000.0]

    chromatograms = []
    chromatograms.append(make_chromatogram(0, "TIC", "TIC", time1, int1))
    chromatograms.append(
        make_chromatogram(
            1,
            "SRM_500.0_250.0",
            "SRM",
            time2,
            int2,
            precursor_mz=500.0,
            product_mz=250.0,
        )
    )
    chromatograms.append(make_chromatogram(2, "BPC", "BPC", time3, int3))

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="0" defaultDataProcessingRef="dp">\n'
    content += "  </spectrumList>\n"
    content += f'  <chromatogramList count="3" defaultDataProcessingRef="dp">\n'
    content += "\n".join(chromatograms) + "\n"
    content += "  </chromatogramList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "with_chromatograms.mzML"), "w") as f:
        f.write(content)


def generate_param_groups():
    """Uses referenceableParamGroupRef for binary encoding params."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [1000.0, 5000.0, 2000.0]

    mz_encoded = encode_doubles_64(mz1)
    int_encoded = encode_doubles_64(int1)

    content = """<?xml version="1.0" encoding="utf-8"?>
<indexedmzML xmlns="http://psi.hupo.org/ms/mzml"
             xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<mzML xmlns="http://psi.hupo.org/ms/mzml">
<cvList count="2">
  <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology"/>
  <cv id="UO" fullName="Unit Ontology"/>
</cvList>
<fileDescription>
  <fileContent>
    <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum"/>
  </fileContent>
</fileDescription>
<referenceableParamGroupList count="2">
  <referenceableParamGroup id="mz_params">
    <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>
    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>
    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>
  </referenceableParamGroup>
  <referenceableParamGroup id="int_params">
    <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>
    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>
    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>
  </referenceableParamGroup>
</referenceableParamGroupList>
"""
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="1" defaultDataProcessingRef="dp">\n'
    content += f'    <spectrum index="0" id="scan=1" defaultArrayLength="3">\n'
    content += '      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>\n'
    content += '      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>\n'
    content += '      <scanList count="1">\n'
    content += "        <scan>\n"
    content += '          <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1.0" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>\n'
    content += "        </scan>\n"
    content += "      </scanList>\n"
    content += '      <binaryDataArrayList count="2">\n'
    content += f'        <binaryDataArray encodedLength="{len(mz_encoded)}">\n'
    content += '          <referenceableParamGroupRef ref="mz_params"/>\n'
    content += f"          <binary>{mz_encoded}</binary>\n"
    content += "        </binaryDataArray>\n"
    content += f'        <binaryDataArray encodedLength="{len(int_encoded)}">\n'
    content += '          <referenceableParamGroupRef ref="int_params"/>\n'
    content += f"          <binary>{int_encoded}</binary>\n"
    content += "        </binaryDataArray>\n"
    content += "      </binaryDataArrayList>\n"
    content += "    </spectrum>\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "param_groups.mzML"), "w") as f:
        f.write(content)


def generate_orphan_ms2():
    """MS2 spectra with no preceding MS1 — ms1_scan_index should be NULL."""
    mz1 = [150.0, 250.0]
    int1 = [3000.0, 1500.0]

    mz2 = [160.0, 260.0]
    int2 = [2000.0, 1000.0]

    spectra = []
    # Spectrum 0: MS2, no preceding MS1
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=2,
            rt_value=1.0,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            mz_values=mz1,
            intensity_values=int1,
            precursor_mz=200.0,
            precursor_charge=2,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    # Spectrum 1: another MS2, still no preceding MS1
    spectra.append(
        make_spectrum(
            index=1,
            spectrum_id="scan=2",
            ms_level=2,
            rt_value=1.1,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            mz_values=mz2,
            intensity_values=int2,
            precursor_mz=300.0,
            precursor_charge=3,
            activation_method="HCD",
            collision_energy=30.0,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="2" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "orphan_ms2.mzML"), "w") as f:
        f.write(content)


def generate_no_index():
    """Bare <mzML> root with no <indexedmzML> wrapper."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [1000.0, 5000.0, 2000.0]

    spectra = []
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.5,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            mz_values=mz1,
            intensity_values=int1,
        )
    )

    # No indexedmzML wrapper — just bare mzML
    content = """<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
<cvList count="2">
  <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology"/>
  <cv id="UO" fullName="Unit Ontology"/>
</cvList>
<fileDescription>
  <fileContent>
    <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum"/>
  </fileContent>
</fileDescription>
"""
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="1" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += "</mzML>\n"

    with open(os.path.join(OUTPUT_DIR, "no_index.mzML"), "w") as f:
        f.write(content)


def generate_mixed_precision():
    """32-bit mz array + 64-bit intensity array in same spectrum."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [1000.0, 5000.0, 2000.0]

    mz_encoded = encode_floats_32(mz1)
    int_encoded = encode_doubles_64(int1)

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="1" defaultDataProcessingRef="dp">\n'
    content += f'    <spectrum index="0" id="scan=1" defaultArrayLength="3">\n'
    content += '      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>\n'
    content += '      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>\n'
    content += '      <scanList count="1">\n'
    content += "        <scan>\n"
    content += '          <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1.0" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>\n'
    content += "        </scan>\n"
    content += "      </scanList>\n"
    content += '      <binaryDataArrayList count="2">\n'
    # mz: 32-bit, no compression
    content += f'        <binaryDataArray encodedLength="{len(mz_encoded)}">\n'
    content += '          <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>\n'
    content += '          <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"/>\n'
    content += '          <cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>\n'
    content += f"          <binary>{mz_encoded}</binary>\n"
    content += "        </binaryDataArray>\n"
    # intensity: 64-bit, no compression
    content += f'        <binaryDataArray encodedLength="{len(int_encoded)}">\n'
    content += '          <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>\n'
    content += '          <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>\n'
    content += '          <cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>\n'
    content += f"          <binary>{int_encoded}</binary>\n"
    content += "        </binaryDataArray>\n"
    content += "      </binaryDataArrayList>\n"
    content += "    </spectrum>\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "mixed_precision.mzML"), "w") as f:
        f.write(content)


def generate_large_spectrum():
    """Single spectrum with 10,000 data points."""
    import math

    n = 10000
    mz_values = [100.0 + i * 0.1 for i in range(n)]
    int_values = [1000.0 * math.sin(i * 0.01) ** 2 for i in range(n)]

    spectra = []
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.0,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=mz_values[int_values.index(max(int_values))],
            base_peak_intensity=max(int_values),
            tic=sum(int_values),
            mz_values=mz_values,
            intensity_values=int_values,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="1" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "large_spectrum.mzML"), "w") as f:
        f.write(content)


def generate_zero_intensity():
    """Spectrum where all intensities are zero (edge case for normalization)."""
    mz1 = [100.0, 200.0, 300.0]
    int1 = [0.0, 0.0, 0.0]

    spectra = []
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.0,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=100.0,
            base_peak_intensity=0.0,
            tic=0.0,
            mz_values=mz1,
            intensity_values=int1,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="1" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "zero_intensity.mzML"), "w") as f:
        f.write(content)


def generate_isotope_pattern():
    """MS1 spectra with known intensity ratios for isotope pattern matching tests.

    Spectrum 0 (MS1): 500:10000, 501:5000, 502:1000, 600:8000
        At ref=500: +1 ratio=5000/10000=0.50, +2 ratio=1000/10000=0.10
        At ref=600: no peaks at 601 or 602
        -> MATCH for target ratios [0.5, 0.1]

    Spectrum 1 (MS1): 500:10000, 501:100, 502:1000
        At ref=500: +1 ratio=100/10000=0.01, +2 ratio=1000/10000=0.10
        -> NO MATCH: +1 ratio 0.01 outside target 0.5 +/- 30% = [0.35, 0.65]
        (Note: +2 ratio matches 0.1 target, but ALL offsets must match simultaneously)

    Spectrum 2 (MS1): 300:6000, 301:3000, 302:600
        At ref=300: +1 ratio=3000/6000=0.50, +2 ratio=600/6000=0.10
        -> MATCH for target ratios [0.5, 0.1]

    Spectrum 3 (MS2): 500:10000, 501:5000 (prec=600)
        MS2, not MS1 — excluded by isotope pattern macro (MS1 only)
    """
    spectra = []

    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.0,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=500.0,
            base_peak_intensity=10000.0,
            tic=24000.0,
            mz_values=[500.0, 501.0, 502.0, 600.0],
            intensity_values=[10000.0, 5000.0, 1000.0, 8000.0],
        )
    )

    spectra.append(
        make_spectrum(
            index=1,
            spectrum_id="scan=2",
            ms_level=1,
            rt_value=2.0,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=500.0,
            base_peak_intensity=10000.0,
            tic=11100.0,
            mz_values=[500.0, 501.0, 502.0],
            intensity_values=[10000.0, 100.0, 1000.0],
        )
    )

    spectra.append(
        make_spectrum(
            index=2,
            spectrum_id="scan=3",
            ms_level=1,
            rt_value=3.0,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=300.0,
            base_peak_intensity=6000.0,
            tic=9600.0,
            mz_values=[300.0, 301.0, 302.0],
            intensity_values=[6000.0, 3000.0, 600.0],
        )
    )

    spectra.append(
        make_spectrum(
            index=3,
            spectrum_id="scan=4",
            ms_level=2,
            rt_value=3.5,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=500.0,
            base_peak_intensity=10000.0,
            tic=15000.0,
            mz_values=[500.0, 501.0],
            intensity_values=[10000.0, 5000.0],
            precursor_mz=600.0,
            precursor_charge=1,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="4" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "isotope_pattern.mzML"), "w") as f:
        f.write(content)


def generate_prec_prod_link():
    """MS2 spectra for testing mzml_x_prec_prod (X from precursor, match product).

    | Idx | Level | precursor_mz | Peaks (mz: intensity)      | Match delta=358.2871? |
    |-----|-------|-------------|---------------------------|----------------------|
    | 0   | MS1   | -           | 100:1000, 200:5000, 500:3000 | N/A                |
    | 1   | MS2   | 500.0       | 150:3000, 250:1500          | No (500-358.29=141.71, no peak) |
    | 2   | MS2   | 500.0       | 141.71:2500, 250:1500       | Yes (peak at 141.71) |
    | 3   | MS2   | 400.0       | 180:1000, 200:2000          | No (400-358.29=41.71, no peak) |
    | 4   | MS2   | 600.0       | 241.71:100, 300:5000        | Yes but low intensity (2%) |

    Spectrum 4: product intensity = 100/5000 = 2% of base peak.
    min_intensity_pct=3 excludes it, min_intensity_pct=1 includes it.
    """
    spectra = []

    # Spectrum 0: MS1
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.0,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=200.0,
            base_peak_intensity=5000.0,
            tic=9000.0,
            mz_values=[100.0, 200.0, 500.0],
            intensity_values=[1000.0, 5000.0, 3000.0],
        )
    )

    # Spectrum 1: MS2, prec=500, no match for delta=358.2871
    spectra.append(
        make_spectrum(
            index=1,
            spectrum_id="scan=2",
            ms_level=2,
            rt_value=1.1,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=150.0,
            base_peak_intensity=3000.0,
            tic=4500.0,
            mz_values=[150.0, 250.0],
            intensity_values=[3000.0, 1500.0],
            precursor_mz=500.0,
            precursor_charge=2,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    # Spectrum 2: MS2, prec=500, match for delta=358.2871 (500-358.2871=141.7129)
    spectra.append(
        make_spectrum(
            index=2,
            spectrum_id="scan=3",
            ms_level=2,
            rt_value=1.2,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=141.71,
            base_peak_intensity=2500.0,
            tic=4000.0,
            mz_values=[141.71, 250.0],
            intensity_values=[2500.0, 1500.0],
            precursor_mz=500.0,
            precursor_charge=2,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    # Spectrum 3: MS2, prec=400, no match for delta=358.2871
    spectra.append(
        make_spectrum(
            index=3,
            spectrum_id="scan=4",
            ms_level=2,
            rt_value=1.3,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=200.0,
            base_peak_intensity=2000.0,
            tic=3000.0,
            mz_values=[180.0, 200.0],
            intensity_values=[1000.0, 2000.0],
            precursor_mz=400.0,
            precursor_charge=2,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    # Spectrum 4: MS2, prec=600, match but low intensity product
    # 600-358.2871=241.7129, peak at 241.71 with intensity 100
    # base_peak = 5000, so i_norm = 100/5000 = 0.02 (2%)
    spectra.append(
        make_spectrum(
            index=4,
            spectrum_id="scan=5",
            ms_level=2,
            rt_value=1.4,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=300.0,
            base_peak_intensity=5000.0,
            tic=5100.0,
            mz_values=[241.71, 300.0],
            intensity_values=[100.0, 5000.0],
            precursor_mz=600.0,
            precursor_charge=2,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="5" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "prec_prod_link.mzML"), "w") as f:
        f.write(content)


def generate_massdefect_prec():
    """MS2 spectra with known precursor mass defects for mzml_x_prec_massdefect.

    | Idx | Level | precursor_mz | Mass defect | Peaks              |
    |-----|-------|-------------|-------------|--------------------|
    | 0   | MS1   | -           | -           | 100:1000, 200:5000 |
    | 1   | MS2   | 200.15      | 0.15        | 150:3000, 250:1500 |
    | 2   | MS2   | 300.92      | 0.92        | 120:500, 220:2500  |
    | 3   | MS2   | 400.50      | 0.50        | 160:1000, 260:2000 |
    """
    spectra = []

    # Spectrum 0: MS1
    spectra.append(
        make_spectrum(
            index=0,
            spectrum_id="scan=1",
            ms_level=1,
            rt_value=1.0,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=200.0,
            base_peak_intensity=5000.0,
            tic=6000.0,
            mz_values=[100.0, 200.0],
            intensity_values=[1000.0, 5000.0],
        )
    )

    # Spectrum 1: MS2, prec=200.15 (defect=0.15)
    spectra.append(
        make_spectrum(
            index=1,
            spectrum_id="scan=2",
            ms_level=2,
            rt_value=1.1,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=150.0,
            base_peak_intensity=3000.0,
            tic=4500.0,
            mz_values=[150.0, 250.0],
            intensity_values=[3000.0, 1500.0],
            precursor_mz=200.15,
            precursor_charge=2,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    # Spectrum 2: MS2, prec=300.92 (defect=0.92)
    spectra.append(
        make_spectrum(
            index=2,
            spectrum_id="scan=3",
            ms_level=2,
            rt_value=1.2,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=220.0,
            base_peak_intensity=2500.0,
            tic=3000.0,
            mz_values=[120.0, 220.0],
            intensity_values=[500.0, 2500.0],
            precursor_mz=300.92,
            precursor_charge=2,
            activation_method="HCD",
            collision_energy=30.0,
        )
    )

    # Spectrum 3: MS2, prec=400.50 (defect=0.50)
    spectra.append(
        make_spectrum(
            index=3,
            spectrum_id="scan=4",
            ms_level=2,
            rt_value=1.3,
            rt_unit="minute",
            spectrum_type="centroid",
            polarity="positive",
            base_peak_mz=260.0,
            base_peak_intensity=2000.0,
            tic=3000.0,
            mz_values=[160.0, 260.0],
            intensity_values=[1000.0, 2000.0],
            precursor_mz=400.50,
            precursor_charge=2,
            activation_method="CID",
            collision_energy=35.0,
        )
    )

    content = INDEXED_MZML_HEADER
    content += RUN_OPEN + "\n"
    content += '  <spectrumList count="4" defaultDataProcessingRef="dp">\n'
    content += "\n".join(spectra) + "\n"
    content += "  </spectrumList>\n"
    content += RUN_CLOSE + "\n"
    content += INDEXED_MZML_FOOTER

    with open(os.path.join(OUTPUT_DIR, "massdefect_prec.mzML"), "w") as f:
        f.write(content)


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    generate_basic_3spectra()
    generate_compressed()
    generate_float32()
    generate_empty_spectra()
    generate_missing_optional()
    generate_with_chromatograms()
    generate_param_groups()
    generate_orphan_ms2()
    generate_no_index()
    generate_mixed_precision()
    generate_large_spectrum()
    generate_zero_intensity()
    generate_isotope_pattern()
    generate_prec_prod_link()
    generate_massdefect_prec()
    print(f"Generated test data in {OUTPUT_DIR}")
