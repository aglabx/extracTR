
import argparse

def run_it(settings):

    ### step 1. Compute aindex for reads

    ### step 2. Find tandem repeats using circular path in de bruijn graph

    ### step 3. Save results to CSV

    ### step 4. Analyze repeat borders

    ### step 5. Enrich repeats variants
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and analyze tandem repeats from raw DNA sequences.")
    parser.add_argument("-1", "--fastq1", help="Input left reads in FASTQ format.")
    parser.add_argument("-2", "--fastq2", help="Input right reads in  in FASTQ format.")
    parser.add_argument("-o", "--output", help="Output prefix for tandem repeats search results.")
    args = parser.parse_args()

    settings = {
        "fastq1": args.fastq1,
        "fastq2": args.fastq2,
        "output": args.output
    }

    print(args)
    run_it(settings)