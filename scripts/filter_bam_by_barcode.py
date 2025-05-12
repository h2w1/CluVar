import pysam

def load_barcodes(barcode_file):
    with open(barcode_file, 'r') as f:
        return set(line.strip() for line in f)

def filter_bam(input_bam, output_bam, barcodes):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)
    
    for read in bam_in:
        try:
            cb = read.get_tag('CB')  # CB:Z: 태그에서 바코드 추출
            if cb in barcodes:
                bam_out.write(read)
        except KeyError:
            # CB 태그가 없는 리드는 건너뜀
            continue
    
    bam_in.close()
    bam_out.close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Filter BAM by cell barcodes.")
    parser.add_argument('-i', '--input_bam', required=True, help="Input BAM file path.")
    parser.add_argument('-o', '--output_bam', required=True, help="Output BAM file path.")
    parser.add_argument('-b', '--barcode_file', required=True, help="File with list of barcodes.")
    
    args = parser.parse_args()

    barcodes = load_barcodes(args.barcode_file)
    filter_bam(args.input_bam, args.output_bam, barcodes)
    print(f"Filtered BAM saved to {args.output_bam}")
