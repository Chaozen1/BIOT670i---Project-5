from Bio import Entrez
import os
import sys
import time

# Always tell NCBI who you are
Entrez.email = "your.email@example.com"

def get_sra_ids_from_samn(samn_ids):
    """
    Fetches SRA run IDs (SRR) for a list of BioSample accession numbers (SAMN).
    """
    sra_ids = []
    
    # Batch search for efficiency
    id_string = " OR ".join(samn_ids)
    
    try:
        # Search the BioSample database first
        handle = Entrez.esearch(db="biosample", term=id_string, retmax=len(samn_ids))
        biosample_record = Entrez.read(handle)
        handle.close()
        
        biosample_ids = biosample_record["IdList"]
        print(f"Found {len(biosample_ids)} BioSample IDs.")

        if not biosample_ids:
            return []

        # Now, link the BioSample IDs to the SRA database to get SRR IDs
        handle = Entrez.elink(dbfrom="biosample", db="sra", id=",".join(biosample_ids))
        link_record = Entrez.read(handle)
        handle.close()

        for link_set in link_record:
            if "LinkSetDb" in link_set and link_set["LinkSetDb"]:
                for link_info in link_set["LinkSetDb"]:
                    if "Link" in link_info:
                        for link in link_info["Link"]:
                            sra_ids.append(link["Id"])

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
    
    return sra_ids

def write_commands_to_file(sra_ids, filename):
    """
    Writes prefetch and fastq-dump commands to a shell script file.
    """
    # Use "w" mode to create or overwrite the file
    with open(filename, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# This script downloads and converts SRA files to FASTQ.\n")
        f.write("# Please make sure the SRA Toolkit is installed and in your PATH.\n\n")
        
        # Write prefetch command
        f.write(f"prefetch {' '.join(sra_ids)}\n")
        
        # Write fastq-dump command
        f.write(f"fastq-dump --split-files --gzip {' '.join(sra_ids)}\n")
    
    print(f"Commands written to {filename}. You can run it from your terminal using 'bash {filename}'")

def main():
    # Generate lists of accession numbers
    gyrB_samn_accessions = [f"SAMN{i}" for i in range(10964863, 10965439)]
    s16S_samn_accessions = [f"SAMN{i}" for i in range(10970131, 10970567)]

    # Fetch SRA IDs for both amplicon types
    print("Fetching SRA accessions for gyrB data...")
    gyrB_sra_ids = get_sra_ids_from_samn(gyrB_samn_accessions)
    print(f"Found {len(gyrB_sra_ids)} gyrB SRA accessions.")
    
    print("Fetching SRA accessions for 16S data...")
    s16S_sra_ids = get_sra_ids_from_samn(s16S_samn_accessions)
    print(f"Found {len(s16S_sra_ids)} 16S SRA accessions.")
    
    # Write the commands to a file
    write_commands_to_file(gyrB_sra_ids, "download_gyrB_data.sh")
    write_commands_to_file(s16S_sra_ids, "download_16S_data.sh")

if __name__ == "__main__":
    main()