import pysam
import argparse
import os

# 명령행 인자를 파싱하기 위한 ArgumentParser를 생성합니다.
parser = argparse.ArgumentParser(description='Retag reads with their cell barcodes and UMIs')

# 필요한 인자들을 추가합니다.
parser.add_argument('-i', '--input', required=True, help="input SAM/BAM file")
parser.add_argument('-o', '--out', required=True, help="output BAM file")
parser.add_argument("--no_umi", required=False, default="False", help="set True if your bam has no umi tag")
parser.add_argument("--umi_tag", required=False, default="UB", help="set if umi tag in output bam should not be UB")
parser.add_argument("--cell_tag", required=False, default="CB", help="set if cell barcode tag should not be CB")
args = parser.parse_args()

# no_umi 인자를 불리언 값으로 변환합니다.
if args.no_umi == "True":
    args.no_umi = True
elif args.no_umi == "False":
    args.no_umi = False

CELL_TAG = args.cell_tag
UMI_TAG = args.umi_tag

# 입력 파일의 확장자를 확인하여 파일 모드를 결정합니다.
input_extension = os.path.splitext(args.input)[1]
if input_extension == '.sam':
    input_mode = 'r'
elif input_extension == '.bam':
    input_mode = 'rb'
else:
    raise ValueError("Input file must be a .sam or .bam file")

# 입력 SAM/BAM 파일을 엽니다.
bam = pysam.AlignmentFile(args.input, input_mode)

# 출력 BAM 파일을 생성합니다.
bamout = pysam.AlignmentFile(args.out, 'wb', template=bam)

# 각 리드를 순회하며 처리합니다.
for read in bam:
    qname = read.query_name  # 리드의 이름을 가져옵니다.
    tokens = qname.split(";")  # 리드 이름을 ';'로 분할합니다.
    if args.no_umi:
        # UMI가 없는 경우, 토큰 수가 2개여야 합니다.
        assert len(tokens) == 2, f"Read name format incorrect: {qname}"
        read.set_tag(CELL_TAG, tokens[-1])  # 셀 바코드 태그를 설정합니다.
    else:
        # UMI가 있는 경우, 토큰 수가 3개여야 합니다.
        assert len(tokens) == 3, f"Read name format incorrect: {qname}"
        read.set_tag(CELL_TAG, tokens[-2])  # 셀 바코드 태그를 설정합니다.
        read.set_tag(UMI_TAG, tokens[-1])   # UMI 태그를 설정합니다.
    # 수정된 리드를 출력 BAM 파일에 작성합니다.
    bamout.write(read)

# 파일을 닫습니다.
bam.close()
bamout.close()
