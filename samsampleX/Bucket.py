import pysam

from samsampleX.Loader import Loader


class Bucket:
    """
    Buffer that holds a certain number of reads
    """

    def __init__(self, out_bam: Loader):
        self.reads: list[pysam.AlignedSegment] = []
        self.out_bam: Loader = out_bam
        self.max_reads: int = 1000

    def __len__(self) -> int:
        """
        Count contained reads
        """
        return len(self.reads)

    def add_read(self, read: pysam.AlignedSegment | list[pysam.AlignedSegment]) -> None:
        """
        Add read to bucket while under maximum count
        """
        if isinstance(read, pysam.AlignedSegment):
            self.reads.append(read)
        elif isinstance(read, list):
            self.reads.extend(read)

        if len(self.reads) > self.max_reads:
            self.write_reads()

    def write_reads(self) -> None:
        """
        Batch write reads when max_reads is reached
        """
        if len(self.reads) > 0:
            for r in self.reads:
                self.out_bam.bam.write(read=r)
            self.reads.clear()
