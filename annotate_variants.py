import sys
import vcf

def annotate_variants(vcf_file, output_file):
    with open(output_file, 'w') as out:
        out.write("Chromosome\tPosition\tRef\tAlt\tGene\tEffect\n")
        
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        for record in vcf_reader:
            gene = record.INFO.get('GENE', ['Unknown'])[0]
            effect = record.INFO.get('EFFECT', ['Unknown'])[0]
            out.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{gene}\t{effect}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 annotate_variants.py <input.vcf> <output.txt>")
        sys.exit(1)

    annotate_variants(sys.argv[1], sys.argv[2])
