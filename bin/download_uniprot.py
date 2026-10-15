import argparse
import datetime
import gzip
import io
import os
import sys
from urllib.request import Request, urlopen

# ftp.uniprot.org itself didn't respond from where this was developed;
# EBI's mirror carries the exact same taxonomic_divisions layout and files.
UNIPROT_FTP_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/"
    "knowledgebase/taxonomic_divisions/")


def log(msg):
    print(f'[{datetime.datetime.now()}] {msg}', file=sys.stderr)


def fetch(url):
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req) as res:
        return res.read()


def download_gzipped(url, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    log(f'Downloading {url}')
    data = fetch(url)
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
        data = gz.read()
    with open(dest + ".part", "wb") as f:
        f.write(data)
    os.replace(dest + ".part", dest)
    log(f'Wrote {dest}')


def main():
    parser = argparse.ArgumentParser(
        description="Download a UniProt taxonomic division flat-file (e.g. "
                    "uniprot_sprot_plants.dat.gz) from UniProt's FTP site. "
                    "This is ONE file shared by every species in that "
                    "division -- unlike STRING/PlantTFDB, there is no "
                    "per-species variant to pick, so riceKB/uniprotToRDF.py "
                    "filters the species it wants out of this file itself.",
        epilog="Example:\n"
               "  %(prog)s --bank sprot --division plants",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bank", default="sprot", choices=["sprot", "trembl"],
                        help="sprot (Swiss-Prot, reviewed) or trembl "
                             "(unreviewed, much larger per division) "
                             "(default: sprot)")
    parser.add_argument("--division", default="plants",
                        help="UniProt taxonomic division name, as used in "
                             "the UniProt FTP filenames (default: plants)")
    parser.add_argument("-o", "--output-dir", default="data/uniprot",
                        help="directory to write the decompressed .dat file "
                             "into (default: data/uniprot)")
    args = parser.parse_args()

    filename = f"uniprot_{args.bank}_{args.division}.dat"
    url = UNIPROT_FTP_BASE + filename + ".gz"
    dest = os.path.join(args.output_dir, filename)
    download_gzipped(url, dest)


if __name__ == "__main__":
    main()
