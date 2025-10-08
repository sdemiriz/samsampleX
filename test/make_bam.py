import pysam
import os

# Utility to make simple BAM files for testing

read_length: int = 100
interval_count: int = 10
reads_per_interval: int = 10

# Header creation
header = pysam.AlignmentHeader.from_dict(
    {
        "HD": {"VN": "1.5"},
        "PG": [],
        "SQ": [{"SN": "chr1", "LN": 248956422}],
        "RG": [],
    }
)

os.makedirs("test/data", exist_ok=True)
out_file = f"test/data/test-{read_length}bp-{interval_count}count.bam"
out_bam = pysam.AlignmentFile(out_file, "wb", header=header)

# Read creation
total_read_count = 0
for interval in range(interval_count):
    interval_start = read_length + interval * read_length + 1

    for i in range(reads_per_interval):
        read = {
            "name": f"TEST_READ_{total_read_count:03d}",
            "ref_name": "chr1",
            "flag": "0",
            "ref_pos": str(interval_start),
            "next_ref_pos": str(interval_start + 1),
            "next_ref_name": "chr1",
            "cigar": f"{read_length}M",
            "seq": "A" * read_length,
            "qual": "?" * read_length,
            "tags": [],
            "map_quality": "60",
            "length": str(read_length * interval_count),
        }

        new_read = pysam.AlignedSegment.from_dict(read, header=header)
        out_bam.write(new_read)

        total_read_count += 1

out_bam.close()
pysam.index(out_file)
