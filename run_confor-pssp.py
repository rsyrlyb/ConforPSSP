from conforpssp import PSSP
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    "--sequence",
    type=str,
    default="L"*50,
    help="target amino acid sequence",
)

parser.add_argument(
    "--fasta_file",
    type=str,
    default=None,
    help="name of a file with target amino acid sequences",
)

parser.add_argument(
    "--output_dir",
    type=str,
    default=None,
    help="directory where to save the output files",
)

parser.add_argument(
    "--N",
    type=int,
    default=15,
    help="number of predicted secondary structures",
)

parser.add_argument(
    "--Y",
    type=int,
    default=15,
    help="how many tokens to consider per position in frequency sequence",
)

parser.add_argument(
    "--model",
    type=int,
    default=1,
    help="model number",
)

args = parser.parse_args()

generator = PSSP(str(args.model))
generator.generate(args.sequence, args.fasta_file, args.output_dir, args.N, args.Y)

